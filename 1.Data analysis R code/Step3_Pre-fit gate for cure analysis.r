# 🔴 Configure: full Step3 pre-fit gate paths, design, and file budget ===============================

## 🟠 Declare: Step1 input folder and Step3 output folder ===============================
# Step1 결과가 저장된 폴더
step1_dir <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦1.Step1_KM/attachments'

# Step3 결과를 저장할 폴더
export_dir <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦3.Step3_non-zero cure fraction screening/attachments'

## 🟠 Declare: study design, dataset targets, and screening thresholds ===============================
# observational: Xie(2024) primary + hdcuremodels::sufficient_fu_test sensitivity
# rct: RECeUS-AIC arm-wise primary
study_design <- 'observational'
# study_design <- 'rct'

# RCT일 때 Step1 bundle$step1_input 안에 arm_var가 보존되어 있어야 한다.
arm_var <- 'arm'

# NULL이면 Step1 summary에 존재하는 모든 dataset prefix를 자동 사용한다.
datasets_to_run <- c('merged', 'pnu', 'snu')
# datasets_to_run <- NULL

alpha_primary <- 0.05
reps_nonzero_cure <- 5000L
reps_xie_boot <- 1000L
min_valid_xie_boot_prop <- 0.80
seed_base <- 20260306L

# NULL이면 hdcuremodels::nonzerocure_test()의 패키지 기본 설정을 사용한다.
b_nonzerocure <- NULL

# RECeUS-AIC 관련 설정
receus_pi_threshold <- 0.025
receus_r_threshold <- 0.05
receus_link <- 'logistic'
receus_candidate_dists <- c('exp', 'weibull', 'gamma', 'llogis')
rct_gate_rule <- 'at_least_one_arm'

plot_width_in <- 8
plot_height_in <- 6
plot_dpi <- 300
max_generated_files <- 20L

# 🔴 Attach: packages and create export folder ===============================

## 🟠 Load: required package set conditional on study design ===============================
required_pkgs <- c(
  'survival', 'hdcuremodels', 'readr', 'dplyr', 'tibble', 'ggplot2'
)

if (identical(study_design, 'rct')) {
  required_pkgs <- c(required_pkgs, 'flexsurv', 'flexsurvcure')
}

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0L) {
  install.packages(to_install)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# 🔴 Define: scalar guards, I/O helpers, and Step1 file discovery ===============================

## 🟠 Implement: scalar coercion, formatting, and CSV utilities ===============================
write_csv_safe <- function(x, path) {
  readr::write_csv(x, path, na = '')
}

assert_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop('Required file not found: ', path)
  }
}

as_scalar_character <- function(x, name) {
  x <- as.character(x)
  if (length(x) < 1L || is.na(x[1]) || x[1] == '') {
    stop(name, ' must be a non-missing scalar character value.')
  }
  x[1]
}

as_scalar_numeric <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) < 1L || !is.finite(x[1])) {
    stop(name, ' must be a finite scalar numeric value.')
  }
  x[1]
}

safe_num_from_list <- function(x, name) {
  if (is.null(x)) return(NA_real_)
  if (!(name %in% names(x))) return(NA_real_)
  val <- x[[name]]
  if (length(val) == 0L) return(NA_real_)
  suppressWarnings(as.numeric(unname(val))[1])
}

safe_first_matching_num <- function(x, candidates) {
  if (is.null(x)) return(NA_real_)
  nm <- intersect(candidates, names(x))
  if (length(nm) < 1L) return(NA_real_)
  val <- x[[nm[1]]]
  if (length(val) == 0L) return(NA_real_)
  suppressWarnings(as.numeric(unname(val))[1])
}

safe_first_matching_chr <- function(x, candidates) {
  if (is.null(x)) return(NA_character_)
  nm <- intersect(candidates, names(x))
  if (length(nm) < 1L) return(NA_character_)
  val <- x[[nm[1]]]
  if (length(val) == 0L) return(NA_character_)
  as.character(unname(val))[1]
}

format_p <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    'NA',
    ifelse(
      x < 10^(-digits),
      paste0('<', formatC(10^(-digits), format = 'f', digits = digits)),
      formatC(x, format = 'f', digits = digits)
    )
  )
}

format_num <- function(x, digits = 3) {
  ifelse(is.na(x), 'NA', formatC(x, format = 'f', digits = digits))
}

prop_ylim_upper <- function(x) {
  vals <- x[is.finite(x)]
  if (length(vals) == 0L) return(1)
  min(1, max(0.10, max(vals) * 1.20 + 0.05))
}

dataset_seed_offset <- function(prefix) {
  as.integer(sum(utf8ToInt(enc2utf8(prefix))))
}

empty_plot <- function(title_txt, subtitle_txt = NULL, body_txt = 'No plottable results') {
  ggplot2::ggplot() +
    ggplot2::annotate('text', x = 0, y = 0, label = body_txt, size = 5) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title_txt, subtitle = subtitle_txt) +
    ggplot2::theme_void(base_size = 12)
}

## 🟠 Implement: Step1 path constructors, dataset resolution, and file budget ===============================
step1_summary_path <- function(step1_dir) {
  file.path(step1_dir, 'step1_dataset_summary.csv')
}

step1_bundle_path <- function(prefix, step1_dir) {
  file.path(step1_dir, paste0(prefix, '_step1_fit_bundle.rds'))
}

load_step1_summary <- function(step1_dir) {
  path <- step1_summary_path(step1_dir)
  assert_file_exists(path)
  
  df <- readr::read_csv(path, show_col_types = FALSE)
  
  required_cols <- c('dataset_label', 'file_prefix')
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      'Step1 summary CSV is missing required columns: ',
      paste(missing_cols, collapse = ', ')
    )
  }
  
  df
}

resolve_datasets_to_run <- function(step1_summary_df, datasets_to_run = NULL) {
  available_prefixes <- unique(as.character(step1_summary_df$file_prefix))
  
  if (is.null(datasets_to_run)) {
    return(available_prefixes)
  }
  
  datasets_to_run <- as.character(datasets_to_run)
  missing_prefixes <- setdiff(datasets_to_run, available_prefixes)
  
  if (length(missing_prefixes) > 0L) {
    stop(
      'Requested dataset prefixes not found in Step1 summary: ',
      paste(missing_prefixes, collapse = ', '),
      '\nAvailable prefixes: ',
      paste(available_prefixes, collapse = ', ')
    )
  }
  
  datasets_to_run
}

make_target_files <- function(prefixes, export_dir) {
  prefixes <- as.character(prefixes)
  
  per_dataset_files <- unlist(lapply(prefixes, function(px) {
    c(
      paste0(px, '_step3_summary.csv'),
      paste0(px, '_step3_gate.png'),
      paste0(px, '_step3_bundle.rds')
    )
  }), use.names = FALSE)
  
  combined_files <- c(
    'step3_dataset_overview.csv',
    'step3_overview_gate.png',
    'step3_pre_fit_gate_all.rds',
    'step3_sessionInfo.txt',
    'step3_manifest.csv'
  )
  
  file.path(export_dir, c(per_dataset_files, combined_files))
}

## 🟠 Implement: Step1 bundle loading, endpoint checks, and carry-forward extraction ===============================
load_step1_bundle <- function(prefix, step1_dir) {
  bundle_path <- step1_bundle_path(prefix, step1_dir)
  assert_file_exists(bundle_path)
  
  bundle <- readRDS(bundle_path)
  
  if (!('fit' %in% names(bundle))) {
    stop('Invalid Step1 bundle for ', prefix, ': missing $fit')
  }
  if (!('km_transition' %in% names(bundle$fit))) {
    stop('Invalid Step1 bundle for ', prefix, ': missing $fit$km_transition')
  }
  if (!inherits(bundle$fit$km_transition, 'survfit')) {
    stop('Invalid Step1 bundle for ', prefix, ': $fit$km_transition is not a survfit object')
  }
  
  bundle
}

check_step1_primary_endpoint <- function(bundle, prefix) {
  if (!('meta' %in% names(bundle))) {
    stop('Invalid Step1 bundle for ', prefix, ': missing $meta')
  }
  if (!('primary_endpoint_definition' %in% names(bundle$meta))) {
    stop('Invalid Step1 bundle for ', prefix, ': missing $meta$primary_endpoint_definition')
  }
  
  ped <- bundle$meta$primary_endpoint_definition
  
  if (!('event' %in% names(ped)) || !('censoring' %in% names(ped))) {
    stop('Invalid Step1 bundle for ', prefix, ': incomplete primary_endpoint_definition')
  }
  
  event_txt <- as.character(ped$event)[1]
  censor_txt <- as.character(ped$censoring)[1]
  
  if (!identical(event_txt, 'transition (status_num == 1)')) {
    stop(
      'Step1 bundle for ', prefix,
      ' is not aligned with the expected primary event definition.\n',
      "Expected event = 'transition (status_num == 1)', found: ", event_txt
    )
  }
  
  if (!identical(censor_txt, 'right_censoring (0) + remission (2)')) {
    stop(
      'Step1 bundle for ', prefix,
      ' is not aligned with the expected censoring definition.\n',
      "Expected censoring = 'right_censoring (0) + remission (2)', found: ", censor_txt
    )
  }
  
  invisible(TRUE)
}

get_step1_summary_row <- function(step1_summary_df, prefix) {
  out <- step1_summary_df |>
    dplyr::filter(file_prefix == !!prefix)
  
  if (nrow(out) != 1L) {
    stop(
      "Step1 summary must contain exactly one row for file_prefix = '",
      prefix, "'. Found: ", nrow(out)
    )
  }
  
  out
}

extract_primary_analysis_df <- function(step1_bundle, arm_var = NULL) {
  if (!('step1_input' %in% names(step1_bundle))) {
    stop('Step1 bundle is missing $step1_input')
  }
  
  dat <- step1_bundle$step1_input
  
  required_cols <- c('days_followup', 'status_num')
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop(
      'Step1 bundle step1_input is missing required columns: ',
      paste(missing_cols, collapse = ', ')
    )
  }
  
  out <- dat |>
    dplyr::mutate(
      time = as.numeric(days_followup),
      event = as.integer(status_num == 1L)
    ) |>
    dplyr::filter(!is.na(time), is.finite(time), time >= 0, !is.na(event)) |>
    dplyr::select(dplyr::any_of(c('site_id', 'id', 'site', arm_var)), time, event, status_num)
  
  if (!is.null(arm_var)) {
    if (!(arm_var %in% names(out))) {
      stop(
        "arm_var = '", arm_var,
        "' is not present in Step1 bundle$step1_input.\n",
        'For study_design = "rct", Step1 must carry this variable forward.'
      )
    }
    out <- out |>
      dplyr::mutate(
        arm_value = as.character(.data[[arm_var]])
      )
  }
  
  out |>
    dplyr::arrange(time)
}

build_step3_carryforward <- function(step1_bundle) {
  list(
    meta = step1_bundle$meta,
    step1_input = step1_bundle$step1_input,
    fit = step1_bundle$fit,
    tables = list(
      summary = step1_bundle$tables$summary,
      timepoints_primary = step1_bundle$tables$timepoints_primary,
      timepoints_all = step1_bundle$tables$timepoints_all
    )
  )
}

# 🔴 Define: Step3A non-zero cure signal functions ===============================

## 🟠 Implement: nonzerocure_test wrapper, flag mapping, and Step3A summaries ===============================
classify_step3a_signal <- function(p_value, test_ok, alpha = 0.05) {
  if (!isTRUE(test_ok)) return('error')
  if (is.na(p_value)) return('review')
  if (p_value < alpha) return('positive_signal')
  'no_positive_signal'
}

map_step3a_gate_flag <- function(step3a_signal_flag) {
  dplyr::case_when(
    step3a_signal_flag == 'positive_signal' ~ 'supports_gate',
    step3a_signal_flag == 'no_positive_signal' ~ 'does_not_support_gate',
    TRUE ~ 'review'
  )
}

make_step3a_note <- function(step3a_signal_flag) {
  dplyr::case_when(
    step3a_signal_flag == 'positive_signal' ~
      'Step3A supports a non-zero cure signal, but Step3B must also support the gate.',
    step3a_signal_flag == 'no_positive_signal' ~
      'Step3A does not show a clear non-zero cure signal.',
    step3a_signal_flag == 'error' ~
      'Step3A failed. Inspect the KM object and Step1 bundle.',
    TRUE ~
      'Step3A requires manual review.'
  )
}

run_nonzerocure_once <- function(km_fit, reps, seed, alpha, b = NULL) {
  if (!inherits(km_fit, 'survfit')) {
    stop('run_nonzerocure_once(): km_fit must be a survfit object.')
  }
  
  test_args <- list(
    object = km_fit,
    reps = reps,
    seed = seed,
    plot = FALSE
  )
  
  if (!is.null(b)) {
    test_args$b <- b
  }
  
  out <- tryCatch(
    {
      res <- do.call(hdcuremodels::nonzerocure_test, test_args)
      list(
        ok = TRUE,
        result = res,
        error_message = NA_character_
      )
    },
    error = function(e) {
      list(
        ok = FALSE,
        result = NULL,
        error_message = conditionMessage(e)
      )
    }
  )
  
  res <- out$result
  p_val <- safe_num_from_list(res, 'p_value')
  signal_flag <- classify_step3a_signal(p_value = p_val, test_ok = isTRUE(out$ok), alpha = alpha)
  gate_flag <- map_step3a_gate_flag(signal_flag)
  
  summary_tbl <- tibble::tibble(
    step3a_test_name = 'hdcuremodels::nonzerocure_test',
    step3a_test_ok = isTRUE(out$ok),
    step3a_proportion_susceptible = safe_num_from_list(res, 'proportion_susceptible'),
    step3a_proportion_cured = safe_num_from_list(res, 'proportion_cured'),
    step3a_p_value = p_val,
    step3a_time_95_percent_of_events_days = safe_num_from_list(res, 'time_95_percent_of_events'),
    step3a_time_95_percent_of_events_years = safe_num_from_list(res, 'time_95_percent_of_events') / 365.25,
    step3a_error_message = out$error_message,
    step3a_signal_flag = signal_flag,
    step3a_gate_flag = gate_flag,
    step3a_note = make_step3a_note(signal_flag)
  )
  
  list(
    summary = summary_tbl,
    raw = out
  )
}

# 🔴 Define: Step3B observational sufficient-follow-up functions ===============================

## 🟠 Implement: KM-derived ingredients for the Xie 2024 statistic ===============================
km_survival_at_time <- function(km_fit, t) {
  t <- as.numeric(t)[1]
  
  if (!is.finite(t)) return(NA_real_)
  if (t <= 0) return(1)
  
  sm <- tryCatch(
    summary(km_fit, times = t, extend = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(sm)) return(NA_real_)
  surv_val <- suppressWarnings(as.numeric(sm$surv)[1])
  
  if (!is.finite(surv_val)) return(NA_real_)
  max(min(surv_val, 1), 0)
}

km_failure_cdf_at_time <- function(km_fit, t) {
  surv_val <- km_survival_at_time(km_fit, t)
  if (!is.finite(surv_val)) return(NA_real_)
  max(min(1 - surv_val, 1), 0)
}

choose_xie_epsilon <- function(t_n, t_k) {
  t_n <- as.numeric(t_n)[1]
  t_k <- as.numeric(t_k)[1]
  
  if (!is.finite(t_n) || !is.finite(t_k)) return(NA_real_)
  if (2 * (t_n - t_k) < t_n) {
    (9 / 8) * t_n - (1 / 4) * t_k
  } else {
    t_n
  }
}

compute_xie_statistic_once <- function(df) {
  if (!all(c('time', 'event') %in% names(df))) {
    stop('compute_xie_statistic_once(): df must contain time and event.')
  }
  
  df <- df |>
    dplyr::filter(!is.na(time), is.finite(time), time >= 0, !is.na(event)) |>
    dplyr::mutate(
      time = as.numeric(time),
      event = as.integer(event)
    )
  
  if (nrow(df) < 2L) {
    return(list(
      ok = FALSE,
      error_message = 'Too few observations for Xie 2024 test.',
      stat = NULL
    ))
  }
  
  if (sum(df$event == 1L, na.rm = TRUE) < 1L) {
    return(list(
      ok = FALSE,
      error_message = 'No observed events. Xie 2024 test is not evaluable.',
      stat = NULL
    ))
  }
  
  km_fit <- tryCatch(
    survival::survfit(survival::Surv(time, event) ~ 1, data = df),
    error = function(e) NULL
  )
  
  if (is.null(km_fit)) {
    return(list(
      ok = FALSE,
      error_message = 'Failed to fit KM object inside Xie 2024 test.',
      stat = NULL
    ))
  }
  
  t_n <- max(df$time, na.rm = TRUE)
  t_k <- max(df$time[df$event == 1L], na.rm = TRUE)
  eps <- choose_xie_epsilon(t_n = t_n, t_k = t_k)
  
  if (!is.finite(eps) || eps <= 0 || eps > t_n) {
    return(list(
      ok = FALSE,
      error_message = 'Invalid epsilon selected for Xie 2024 test.',
      stat = NULL
    ))
  }
  
  p_hat_n <- km_failure_cdf_at_time(km_fit, t_n)
  f_eps <- km_failure_cdf_at_time(km_fit, t_n - eps)
  f_half <- km_failure_cdf_at_time(km_fit, t_n - eps / 2)
  
  denom <- 2 * f_half - f_eps - p_hat_n
  
  if (!is.finite(denom) || abs(denom) < .Machine$double.eps) {
    p_hat_g <- p_hat_n
  } else {
    p_hat_g <- f_eps + (f_half - f_eps)^2 / denom
  }
  
  if (!is.finite(p_hat_g) || p_hat_g < p_hat_n) {
    p_hat_g <- p_hat_n
  }
  if (p_hat_g > 1) {
    p_hat_g <- 1
  }
  
  T_n <- p_hat_g - p_hat_n
  if (!is.finite(T_n) || T_n < 0) {
    T_n <- 0
  }
  
  list(
    ok = TRUE,
    error_message = NA_character_,
    stat = list(
      n = nrow(df),
      t_n = t_n,
      t_k = t_k,
      epsilon = eps,
      p_hat_n = p_hat_n,
      p_hat_g = p_hat_g,
      T_n = T_n
    )
  )
}

## 🟠 Implement: bootstrap-based Xie 2024 runner under H0 sufficient follow-up ===============================
run_xie_test <- function(df,
                         alpha = 0.05,
                         B = 1000L,
                         seed = 1L,
                         min_valid_prop = 0.80) {
  base_stat <- compute_xie_statistic_once(df)
  
  if (!isTRUE(base_stat$ok)) {
    summary_tbl <- tibble::tibble(
      xie_test_name = 'Xie_2024_bootstrap',
      xie_test_ok = FALSE,
      xie_largest_observed_time_days = NA_real_,
      xie_largest_event_time_days = NA_real_,
      xie_epsilon_days = NA_real_,
      xie_p_hat_n = NA_real_,
      xie_p_hat_g = NA_real_,
      xie_Tn = NA_real_,
      xie_bootstrap_B = as.integer(B),
      xie_bootstrap_valid_B = NA_real_,
      xie_critical_value = NA_real_,
      xie_bootstrap_p_value_approx = NA_real_,
      xie_reject_h0_sufficient = NA,
      xie_primary_flag = 'review',
      xie_error_message = base_stat$error_message
    )
    
    return(list(
      summary = summary_tbl,
      raw = list(base = base_stat, T_boot = NULL)
    ))
  }
  
  set.seed(seed)
  
  n <- nrow(df)
  T_n <- base_stat$stat$T_n
  T_boot <- rep(NA_real_, B)
  
  for (b in seq_len(B)) {
    idx <- sample.int(n = n, size = n, replace = TRUE)
    boot_df <- df[idx, , drop = FALSE]
    boot_stat <- compute_xie_statistic_once(boot_df)
    
    if (isTRUE(boot_stat$ok)) {
      T_boot[b] <- boot_stat$stat$T_n
    }
  }
  
  valid_B <- sum(is.finite(T_boot))
  min_valid_B <- ceiling(B * min_valid_prop)
  
  if (valid_B < min_valid_B) {
    summary_tbl <- tibble::tibble(
      xie_test_name = 'Xie_2024_bootstrap',
      xie_test_ok = FALSE,
      xie_largest_observed_time_days = base_stat$stat$t_n,
      xie_largest_event_time_days = base_stat$stat$t_k,
      xie_epsilon_days = base_stat$stat$epsilon,
      xie_p_hat_n = base_stat$stat$p_hat_n,
      xie_p_hat_g = base_stat$stat$p_hat_g,
      xie_Tn = base_stat$stat$T_n,
      xie_bootstrap_B = as.integer(B),
      xie_bootstrap_valid_B = valid_B,
      xie_critical_value = NA_real_,
      xie_bootstrap_p_value_approx = NA_real_,
      xie_reject_h0_sufficient = NA,
      xie_primary_flag = 'review',
      xie_error_message = paste0(
        'Too few valid bootstrap replicates for Xie 2024 test: ',
        valid_B, ' / ', B
      )
    )
    
    return(list(
      summary = summary_tbl,
      raw = list(base = base_stat, T_boot = T_boot)
    ))
  }
  
  centered_boot <- T_boot - T_n
  crit_val <- as.numeric(stats::quantile(centered_boot, probs = 1 - alpha, na.rm = TRUE, type = 8))
  reject_h0 <- isTRUE(T_n > crit_val)
  p_val_approx <- mean(centered_boot >= T_n, na.rm = TRUE)
  
  primary_flag <- dplyr::case_when(
    reject_h0 ~ 'does_not_support_gate',
    !reject_h0 ~ 'supports_gate',
    TRUE ~ 'review'
  )
  
  summary_tbl <- tibble::tibble(
    xie_test_name = 'Xie_2024_bootstrap',
    xie_test_ok = TRUE,
    xie_largest_observed_time_days = base_stat$stat$t_n,
    xie_largest_event_time_days = base_stat$stat$t_k,
    xie_epsilon_days = base_stat$stat$epsilon,
    xie_p_hat_n = base_stat$stat$p_hat_n,
    xie_p_hat_g = base_stat$stat$p_hat_g,
    xie_Tn = base_stat$stat$T_n,
    xie_bootstrap_B = as.integer(B),
    xie_bootstrap_valid_B = valid_B,
    xie_critical_value = crit_val,
    xie_bootstrap_p_value_approx = p_val_approx,
    xie_reject_h0_sufficient = reject_h0,
    xie_primary_flag = primary_flag,
    xie_error_message = NA_character_
  )
  
  list(
    summary = summary_tbl,
    raw = list(base = base_stat, T_boot = T_boot)
  )
}

## 🟠 Implement: hdcuremodels sufficient follow-up sensitivity wrapper ===============================
run_mz_sufficient_fu_test <- function(km_fit, alpha = 0.05) {
  if (!inherits(km_fit, 'survfit')) {
    stop('run_mz_sufficient_fu_test(): km_fit must be a survfit object.')
  }
  
  out <- tryCatch(
    {
      res <- hdcuremodels::sufficient_fu_test(object = km_fit)
      list(
        ok = TRUE,
        result = res,
        error_message = NA_character_
      )
    },
    error = function(e) {
      list(
        ok = FALSE,
        result = NULL,
        error_message = conditionMessage(e)
      )
    }
  )
  
  res <- out$result
  p_val <- safe_first_matching_num(res, c('p_value', 'p.value'))
  stat_val <- safe_first_matching_num(res, c('q_n', 'Q_n', 'test_statistic', 'statistic'))
  
  sensitivity_flag <- dplyr::case_when(
    !isTRUE(out$ok) ~ 'review',
    is.na(p_val) ~ 'review',
    p_val < alpha ~ 'supports_gate',
    TRUE ~ 'does_not_support_gate'
  )
  
  summary_tbl <- tibble::tibble(
    mz_test_name = 'hdcuremodels::sufficient_fu_test',
    mz_test_ok = isTRUE(out$ok),
    mz_statistic = stat_val,
    mz_p_value = p_val,
    mz_supports_sufficient_followup = ifelse(is.na(p_val), NA, p_val < alpha),
    mz_sensitivity_flag = sensitivity_flag,
    mz_error_message = out$error_message
  )
  
  list(
    summary = summary_tbl,
    raw = out
  )
}

make_step3b_observational_note <- function(primary_flag, sensitivity_flag) {
  if (identical(primary_flag, 'supports_gate') && identical(sensitivity_flag, 'supports_gate')) {
    return('Step3B observational primary and sensitivity both support sufficient follow-up.')
  }
  if (identical(primary_flag, 'supports_gate') && identical(sensitivity_flag, 'does_not_support_gate')) {
    return('Step3B observational primary supports sufficient follow-up, but the sensitivity check does not.')
  }
  if (identical(primary_flag, 'does_not_support_gate')) {
    return('Step3B observational primary does not support sufficient follow-up.')
  }
  'Step3B observational results require manual review.'
}

run_step3b_observational <- function(df,
                                     km_fit,
                                     alpha,
                                     reps_xie_boot,
                                     seed_xie,
                                     min_valid_prop) {
  xie_res <- run_xie_test(
    df = df,
    alpha = alpha,
    B = reps_xie_boot,
    seed = seed_xie,
    min_valid_prop = min_valid_prop
  )
  
  mz_res <- run_mz_sufficient_fu_test(
    km_fit = km_fit,
    alpha = alpha
  )
  
  primary_flag <- xie_res$summary$xie_primary_flag[1]
  sensitivity_flag <- mz_res$summary$mz_sensitivity_flag[1]
  
  summary_tbl <- tibble::tibble(
    step3b_primary_method = 'Xie_2024_bootstrap',
    step3b_sensitivity_method = 'hdcuremodels::sufficient_fu_test',
    step3b_primary_gate_flag = primary_flag,
    step3b_sensitivity_gate_flag = sensitivity_flag,
    step3b_note = make_step3b_observational_note(primary_flag, sensitivity_flag)
  ) |>
    dplyr::bind_cols(xie_res$summary) |>
    dplyr::bind_cols(mz_res$summary)
  
  list(
    summary = summary_tbl,
    raw = list(
      xie = xie_res$raw,
      mz = mz_res$raw
    )
  )
}

# 🔴 Define: Step3B RCT RECeUS-AIC functions ===============================

## 🟠 Implement: safe parametric model fitting helpers for RECeUS-AIC ===============================
extract_aic_from_fit <- function(fit) {
  val <- tryCatch({
    if (!is.null(fit$AIC)) {
      as.numeric(fit$AIC)
    } else {
      as.numeric(stats::AIC(fit))
    }
  }, error = function(e) NA_real_)
  
  if (!is.finite(val)) return(NA_real_)
  val[1]
}

extract_survival_at_time <- function(fit, t) {
  sm <- tryCatch(
    summary(fit, t = t, type = 'survival', ci = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(sm)) return(NA_real_)
  
  if (is.data.frame(sm) && 'est' %in% names(sm)) {
    return(as.numeric(sm$est[1]))
  }
  
  if (is.list(sm) && length(sm) >= 1L) {
    first <- sm[[1]]
    if (is.data.frame(first) && 'est' %in% names(first)) {
      return(as.numeric(first$est[1]))
    }
  }
  
  NA_real_
}

extract_cure_fraction_logistic <- function(fit) {
  coef_vec <- tryCatch(stats::coef(fit), error = function(e) NULL)
  if (is.null(coef_vec)) return(NA_real_)
  
  nm <- names(coef_vec)
  idx <- grep('^theta($|[^[:alnum:]_])|^theta|^pi($|[^[:alnum:]_])|^pi_', nm)
  if (length(idx) < 1L) {
    idx <- grep('theta|pi', nm)
  }
  if (length(idx) < 1L) return(NA_real_)
  
  eta <- suppressWarnings(as.numeric(coef_vec[idx[1]]))
  if (!is.finite(eta)) return(NA_real_)
  
  plogis(eta)
}

safe_fit_noncure_model <- function(df, dist_name) {
  out <- tryCatch(
    {
      fit <- flexsurv::flexsurvreg(
        formula = survival::Surv(time, event) ~ 1,
        data = df,
        dist = dist_name
      )
      
      list(
        ok = TRUE,
        fit = fit,
        error_message = NA_character_
      )
    },
    error = function(e) {
      list(
        ok = FALSE,
        fit = NULL,
        error_message = conditionMessage(e)
      )
    }
  )
  
  list(
    ok = isTRUE(out$ok),
    fit = out$fit,
    error_message = out$error_message,
    model_class = 'noncure',
    model_name = paste0('noncure_', dist_name),
    dist = dist_name,
    aic = if (isTRUE(out$ok)) extract_aic_from_fit(out$fit) else NA_real_
  )
}

safe_fit_mixture_cure_model <- function(df, dist_name, link_name = 'logistic') {
  out <- tryCatch(
    {
      fit <- flexsurvcure::flexsurvcure(
        formula = survival::Surv(time, event) ~ 1,
        data = df,
        dist = dist_name,
        link = link_name,
        mixture = TRUE
      )
      
      list(
        ok = TRUE,
        fit = fit,
        error_message = NA_character_
      )
    },
    error = function(e) {
      list(
        ok = FALSE,
        fit = NULL,
        error_message = conditionMessage(e)
      )
    }
  )
  
  list(
    ok = isTRUE(out$ok),
    fit = out$fit,
    error_message = out$error_message,
    model_class = 'mixture_cure',
    model_name = paste0('mixture_cure_', dist_name),
    dist = dist_name,
    aic = if (isTRUE(out$ok)) extract_aic_from_fit(out$fit) else NA_real_
  )
}

## 🟠 Implement: arm-wise RECeUS-AIC model selection, ratio estimation, and thresholds ===============================
fit_receus_aic_one_arm <- function(df_arm,
                                   pi_threshold = 0.025,
                                   r_threshold = 0.05,
                                   candidate_dists = c('exp', 'weibull', 'gamma', 'llogis'),
                                   link_name = 'logistic') {
  if (nrow(df_arm) < 5L) {
    return(list(
      ok = FALSE,
      error_message = 'Too few observations for RECeUS-AIC in this arm.',
      arm_table = tibble::tibble(),
      best_fit = NULL,
      arm_summary = tibble::tibble(
        arm_value = unique(as.character(df_arm$arm_value))[1],
        arm_n = nrow(df_arm),
        arm_events = sum(df_arm$event == 1L, na.rm = TRUE),
        arm_tau_days = max(df_arm$time, na.rm = TRUE),
        arm_best_model_class = NA_character_,
        arm_best_model_name = NA_character_,
        arm_best_aic = NA_real_,
        arm_pi_hat = NA_real_,
        arm_s_overall_tau = NA_real_,
        arm_s_uc_tau = NA_real_,
        arm_r_hat = NA_real_,
        arm_pass = NA,
        arm_note = 'Too few observations for RECeUS-AIC.'
      )
    ))
  }
  
  fit_list <- list()
  
  for (dist_name in candidate_dists) {
    fit_list[[paste0('noncure_', dist_name)]] <- safe_fit_noncure_model(df_arm, dist_name)
    fit_list[[paste0('mixture_cure_', dist_name)]] <- safe_fit_mixture_cure_model(df_arm, dist_name, link_name = link_name)
  }
  
  arm_table <- dplyr::bind_rows(lapply(fit_list, function(x) {
    tibble::tibble(
      model_class = x$model_class,
      model_name = x$model_name,
      dist = x$dist,
      fit_ok = x$ok,
      aic = x$aic,
      error_message = x$error_message
    )
  })) |>
    dplyr::arrange(aic)
  
  valid_names <- names(fit_list)[vapply(fit_list, function(x) isTRUE(x$ok) && is.finite(x$aic), logical(1))]
  
  if (length(valid_names) < 1L) {
    return(list(
      ok = FALSE,
      error_message = 'All RECeUS-AIC candidate fits failed.',
      arm_table = arm_table,
      best_fit = NULL,
      arm_summary = tibble::tibble(
        arm_value = unique(as.character(df_arm$arm_value))[1],
        arm_n = nrow(df_arm),
        arm_events = sum(df_arm$event == 1L, na.rm = TRUE),
        arm_tau_days = max(df_arm$time, na.rm = TRUE),
        arm_best_model_class = NA_character_,
        arm_best_model_name = NA_character_,
        arm_best_aic = NA_real_,
        arm_pi_hat = NA_real_,
        arm_s_overall_tau = NA_real_,
        arm_s_uc_tau = NA_real_,
        arm_r_hat = NA_real_,
        arm_pass = NA,
        arm_note = 'All RECeUS-AIC candidate fits failed.'
      )
    ))
  }
  
  best_name <- valid_names[which.min(vapply(fit_list[valid_names], function(x) x$aic, numeric(1)))]
  best_fit_obj <- fit_list[[best_name]]
  
  tau <- max(df_arm$time, na.rm = TRUE)
  
  if (identical(best_fit_obj$model_class, 'noncure')) {
    arm_summary <- tibble::tibble(
      arm_value = unique(as.character(df_arm$arm_value))[1],
      arm_n = nrow(df_arm),
      arm_events = sum(df_arm$event == 1L, na.rm = TRUE),
      arm_tau_days = tau,
      arm_best_model_class = best_fit_obj$model_class,
      arm_best_model_name = best_fit_obj$model_name,
      arm_best_aic = best_fit_obj$aic,
      arm_pi_hat = 0,
      arm_s_overall_tau = extract_survival_at_time(best_fit_obj$fit, tau),
      arm_s_uc_tau = NA_real_,
      arm_r_hat = NA_real_,
      arm_pass = FALSE,
      arm_note = 'Best AIC model is noncure, so RECeUS-AIC does not support cure-model appropriateness in this arm.'
    )
    
    return(list(
      ok = TRUE,
      error_message = NA_character_,
      arm_table = arm_table,
      best_fit = best_fit_obj,
      arm_summary = arm_summary
    ))
  }
  
  pi_hat <- extract_cure_fraction_logistic(best_fit_obj$fit)
  s_overall_tau <- extract_survival_at_time(best_fit_obj$fit, tau)
  
  if (!is.finite(pi_hat) || !is.finite(s_overall_tau)) {
    arm_summary <- tibble::tibble(
      arm_value = unique(as.character(df_arm$arm_value))[1],
      arm_n = nrow(df_arm),
      arm_events = sum(df_arm$event == 1L, na.rm = TRUE),
      arm_tau_days = tau,
      arm_best_model_class = best_fit_obj$model_class,
      arm_best_model_name = best_fit_obj$model_name,
      arm_best_aic = best_fit_obj$aic,
      arm_pi_hat = pi_hat,
      arm_s_overall_tau = s_overall_tau,
      arm_s_uc_tau = NA_real_,
      arm_r_hat = NA_real_,
      arm_pass = NA,
      arm_note = 'Best mixture-cure model fitted, but pi_hat or S(tau) could not be extracted.'
    )
    
    return(list(
      ok = FALSE,
      error_message = 'Failed to extract pi_hat or S(tau) from best RECeUS-AIC model.',
      arm_table = arm_table,
      best_fit = best_fit_obj,
      arm_summary = arm_summary
    ))
  }
  
  if ((1 - pi_hat) <= .Machine$double.eps) {
    s_uc_tau <- 0
  } else {
    s_uc_tau <- (s_overall_tau - pi_hat) / (1 - pi_hat)
  }
  
  s_uc_tau <- max(min(s_uc_tau, 1), 0)
  r_hat <- ifelse(s_overall_tau > 0, s_uc_tau / s_overall_tau, NA_real_)
  if (is.finite(r_hat)) {
    r_hat <- max(min(r_hat, 1), 0)
  }
  
  arm_pass <- isTRUE(pi_hat > pi_threshold && r_hat < r_threshold)
  
  arm_summary <- tibble::tibble(
    arm_value = unique(as.character(df_arm$arm_value))[1],
    arm_n = nrow(df_arm),
    arm_events = sum(df_arm$event == 1L, na.rm = TRUE),
    arm_tau_days = tau,
    arm_best_model_class = best_fit_obj$model_class,
    arm_best_model_name = best_fit_obj$model_name,
    arm_best_aic = best_fit_obj$aic,
    arm_pi_hat = pi_hat,
    arm_s_overall_tau = s_overall_tau,
    arm_s_uc_tau = s_uc_tau,
    arm_r_hat = r_hat,
    arm_pass = arm_pass,
    arm_note = ifelse(
      arm_pass,
      'This arm passes RECeUS-AIC thresholds.',
      'This arm does not pass RECeUS-AIC thresholds.'
    )
  )
  
  list(
    ok = TRUE,
    error_message = NA_character_,
    arm_table = arm_table,
    best_fit = best_fit_obj,
    arm_summary = arm_summary
  )
}

make_step3b_rct_note <- function(primary_flag, pass_arms) {
  if (identical(primary_flag, 'supports_gate')) {
    paste0('At least one arm passes RECeUS-AIC. Passing arms: ', pass_arms)
  } else if (identical(primary_flag, 'does_not_support_gate')) {
    'No arm passes RECeUS-AIC.'
  } else {
    'Step3B RCT results require manual review.'
  }
}

run_step3b_rct <- function(df,
                           arm_var,
                           pi_threshold = 0.025,
                           r_threshold = 0.05,
                           candidate_dists = c('exp', 'weibull', 'gamma', 'llogis'),
                           link_name = 'logistic',
                           gate_rule = 'at_least_one_arm') {
  if (!(arm_var %in% names(df))) {
    stop(
      "arm_var = '", arm_var,
      "' is not available in the Step1 carry-forward data.\n",
      'For study_design = "rct", Step1 must preserve this variable.'
    )
  }
  
  df <- df |>
    dplyr::mutate(
      arm_value = as.character(.data[[arm_var]])
    )
  
  arm_split <- split(df, df$arm_value)
  arm_results <- lapply(arm_split, function(d) {
    fit_receus_aic_one_arm(
      df_arm = d,
      pi_threshold = pi_threshold,
      r_threshold = r_threshold,
      candidate_dists = candidate_dists,
      link_name = link_name
    )
  })
  
  arm_summary_tbl <- dplyr::bind_rows(lapply(arm_results, `[[`, 'arm_summary'))
  
  n_arms_total <- length(arm_split)
  n_arms_evaluable <- sum(!is.na(arm_summary_tbl$arm_pass))
  any_arm_pass <- any(arm_summary_tbl$arm_pass %in% TRUE, na.rm = TRUE)
  pass_arms <- arm_summary_tbl$arm_value[arm_summary_tbl$arm_pass %in% TRUE]
  
  primary_flag <- dplyr::case_when(
    n_arms_evaluable < 1L ~ 'review',
    identical(gate_rule, 'at_least_one_arm') && any_arm_pass ~ 'supports_gate',
    identical(gate_rule, 'at_least_one_arm') && !any_arm_pass ~ 'does_not_support_gate',
    TRUE ~ 'review'
  )
  
  best_reporting_row <- NULL
  
  if (any_arm_pass) {
    best_reporting_row <- arm_summary_tbl |>
      dplyr::filter(arm_pass %in% TRUE) |>
      dplyr::arrange(arm_r_hat, arm_best_aic) |>
      dplyr::slice(1)
  } else {
    best_reporting_row <- arm_summary_tbl |>
      dplyr::arrange(arm_best_aic) |>
      dplyr::slice(1)
  }
  
  summary_tbl <- tibble::tibble(
    step3b_primary_method = 'RECeUS_AIC_armwise',
    step3b_sensitivity_method = NA_character_,
    step3b_primary_gate_flag = primary_flag,
    step3b_sensitivity_gate_flag = NA_character_,
    step3b_note = make_step3b_rct_note(
      primary_flag = primary_flag,
      pass_arms = ifelse(length(pass_arms) < 1L, 'none', paste(pass_arms, collapse = '; '))
    ),
    receus_gate_rule = gate_rule,
    receus_pi_threshold = pi_threshold,
    receus_r_threshold = r_threshold,
    receus_arm_var = arm_var,
    receus_n_arms_total = n_arms_total,
    receus_n_arms_evaluable = n_arms_evaluable,
    receus_any_arm_pass = any_arm_pass,
    receus_pass_arms = ifelse(length(pass_arms) < 1L, 'none', paste(pass_arms, collapse = '; ')),
    receus_best_arm = best_reporting_row$arm_value[1],
    receus_best_model_class = best_reporting_row$arm_best_model_class[1],
    receus_best_model_name = best_reporting_row$arm_best_model_name[1],
    receus_best_aic = best_reporting_row$arm_best_aic[1],
    receus_best_pi_hat = best_reporting_row$arm_pi_hat[1],
    receus_best_s_overall_tau = best_reporting_row$arm_s_overall_tau[1],
    receus_best_s_uc_tau = best_reporting_row$arm_s_uc_tau[1],
    receus_best_r_hat = best_reporting_row$arm_r_hat[1],
    receus_best_tau_days = best_reporting_row$arm_tau_days[1],
    receus_error_message = NA_character_
  )
  
  list(
    summary = summary_tbl,
    raw = list(
      arm_results = arm_results,
      arm_summary = arm_summary_tbl
    )
  )
}

# 🔴 Define: gate classification, summary construction, and CSV-driven plots ===============================

## 🟠 Implement: overall Step3 gate classification and summary assembly ===============================
classify_step3_gate <- function(step3a_gate_flag,
                                step3b_primary_gate_flag,
                                step3b_sensitivity_gate_flag = NA_character_) {
  if (any(c(step3a_gate_flag, step3b_primary_gate_flag) %in% 'review')) {
    return('review')
  }
  
  a_support <- identical(step3a_gate_flag, 'supports_gate')
  b_support <- identical(step3b_primary_gate_flag, 'supports_gate')
  
  sens_available <- !is.na(step3b_sensitivity_gate_flag)
  sens_support <- identical(step3b_sensitivity_gate_flag, 'supports_gate')
  
  if (a_support && b_support && (!sens_available || sens_support)) {
    return('pass')
  }
  
  if ((a_support && b_support && sens_available && !sens_support) ||
      xor(a_support, b_support)) {
    return('borderline')
  }
  
  'fail'
}

make_step3_gate_note <- function(step3_gate) {
  dplyr::case_when(
    step3_gate == 'pass' ~
      'Step3 gate passes. Cure-model fitting can proceed to the next step.',
    step3_gate == 'borderline' ~
      'Step3 gate is borderline. Cure models should be treated as sensitivity analysis unless later evidence is strong.',
    step3_gate == 'fail' ~
      'Step3 gate fails. Cure models should not be used as the primary analysis at this stage.',
    TRUE ~
      'Step3 gate requires manual review.'
  )
}

build_step3_summary <- function(step1_row,
                                study_design,
                                step3a_tbl,
                                step3b_tbl,
                                alpha,
                                reps_nonzero_cure,
                                reps_xie_boot,
                                seed_nonzero_cure,
                                seed_xie) {
  step1_keep <- step1_row |>
    dplyr::select(dplyr::any_of(c(
      'dataset_label',
      'file_prefix',
      'n',
      'n_right_censoring',
      'n_transition',
      'n_remission',
      'n_primary_events',
      'n_primary_censored',
      'min_followup_days',
      'median_observed_time_days_raw',
      'max_followup_days',
      'max_followup_years',
      'reversekm_median_days',
      'reversekm_median_lcl95_days',
      'reversekm_median_ucl95_days',
      'reversekm_rmean_days',
      'reversekm_se_rmean_days',
      'last_transition_days',
      'plateau_length_days',
      'n_risk_at_last_transition',
      'km_at_last_transition',
      'one_more_event_drop_pct',
      'n_right_censor_after_last_transition',
      'n_remission_after_last_transition',
      'n_primary_censor_after_last_transition',
      'one_minus_km_interpretation',
      'reverse_km_interpretation',
      'reporting_rule',
      'betensky_stability_limits_implemented'
    )))
  
  step3_gate <- classify_step3_gate(
    step3a_gate_flag = step3a_tbl$step3a_gate_flag[1],
    step3b_primary_gate_flag = step3b_tbl$step3b_primary_gate_flag[1],
    step3b_sensitivity_gate_flag = step3b_tbl$step3b_sensitivity_gate_flag[1]
  )
  
  tibble::tibble(
    analysis_step = 'Step3_full_pre_fit_gate',
    study_design = study_design,
    alpha = alpha,
    reps_nonzero_cure = reps_nonzero_cure,
    reps_xie_boot = if (identical(study_design, 'observational')) reps_xie_boot else NA_real_,
    seed_nonzero_cure = seed_nonzero_cure,
    seed_xie = if (identical(study_design, 'observational')) seed_xie else NA_real_
  ) |>
    dplyr::bind_cols(step1_keep) |>
    dplyr::bind_cols(step3a_tbl) |>
    dplyr::bind_cols(step3b_tbl) |>
    dplyr::mutate(
      step3_gate = step3_gate,
      step3_gate_note = make_step3_gate_note(step3_gate)
    )
}

## 🟠 Implement: per-dataset and overview gate plots from exported CSV files ===============================
render_step3_dataset_gate_plot_from_csv <- function(summary_csv, out_png, title_txt) {
  df <- readr::read_csv(summary_csv, show_col_types = FALSE)
  
  if (nrow(df) == 0L) {
    p <- empty_plot(
      title_txt = title_txt,
      subtitle_txt = 'Gate summary',
      body_txt = 'No Step3 summary row available'
    )
  } else {
    sensitivity_flag <- df$step3b_sensitivity_gate_flag[1]
    if (is.na(sensitivity_flag) || sensitivity_flag == '') {
      sensitivity_flag <- 'not_applicable'
    }
    
    plot_df <- tibble::tibble(
      stage = c('Step3A', 'Step3B_primary', 'Step3B_sensitivity', 'Overall_gate'),
      flag = c(
        df$step3a_gate_flag[1],
        df$step3b_primary_gate_flag[1],
        sensitivity_flag,
        df$step3_gate[1]
      )
    ) |>
      dplyr::mutate(
        stage = factor(stage, levels = c('Step3A', 'Step3B_primary', 'Step3B_sensitivity', 'Overall_gate'))
      )
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = stage, y = 1, fill = flag)) +
      ggplot2::geom_tile(height = 0.7, width = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = flag), size = 3.4) +
      ggplot2::scale_fill_manual(
        values = c(
          supports_gate = '#1b9e77',
          does_not_support_gate = '#d95f02',
          review = '#7570b3',
          not_applicable = '#bdbdbd',
          pass = '#1b9e77',
          borderline = '#e6ab02',
          fail = '#d95f02'
        ),
        drop = FALSE
      ) +
      ggplot2::scale_y_continuous(NULL, breaks = NULL) +
      ggplot2::labs(
        title = title_txt,
        subtitle = paste0(
          'Design = ', df$study_design[1],
          ' | Step3 gate = ', df$step3_gate[1],
          ' | n = ', df$n[1],
          ' | events = ', df$n_transition[1]
        ),
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = 'none',
        panel.grid = ggplot2::element_blank()
      )
  }
  
  ggplot2::ggsave(
    filename = out_png,
    plot = p,
    width = plot_width_in,
    height = plot_height_in,
    dpi = plot_dpi
  )
}

render_step3_overview_gate_plot_from_csv <- function(overview_csv, out_png, title_txt, dataset_order) {
  df <- readr::read_csv(overview_csv, show_col_types = FALSE)
  
  if (nrow(df) == 0L) {
    p <- empty_plot(
      title_txt = title_txt,
      subtitle_txt = 'Overview gate summary',
      body_txt = 'No Step3 overview rows available'
    )
  } else {
    plot_df <- dplyr::bind_rows(
      df |>
        dplyr::transmute(file_prefix, stage = 'Step3A', flag = step3a_gate_flag),
      df |>
        dplyr::transmute(file_prefix, stage = 'Step3B_primary', flag = step3b_primary_gate_flag),
      df |>
        dplyr::transmute(
          file_prefix,
          stage = 'Step3B_sensitivity',
          flag = ifelse(is.na(step3b_sensitivity_gate_flag) | step3b_sensitivity_gate_flag == '', 'not_applicable', step3b_sensitivity_gate_flag)
        ),
      df |>
        dplyr::transmute(file_prefix, stage = 'Overall_gate', flag = step3_gate)
    ) |>
      dplyr::mutate(
        file_prefix = factor(file_prefix, levels = dataset_order),
        stage = factor(stage, levels = c('Step3A', 'Step3B_primary', 'Step3B_sensitivity', 'Overall_gate'))
      )
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = stage, y = file_prefix, fill = flag)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = flag), size = 3.0) +
      ggplot2::scale_fill_manual(
        values = c(
          supports_gate = '#1b9e77',
          does_not_support_gate = '#d95f02',
          review = '#7570b3',
          not_applicable = '#bdbdbd',
          pass = '#1b9e77',
          borderline = '#e6ab02',
          fail = '#d95f02'
        ),
        drop = FALSE
      ) +
      ggplot2::labs(
        title = title_txt,
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = 'none',
        panel.grid = ggplot2::element_blank()
      )
  }
  
  ggplot2::ggsave(
    filename = out_png,
    plot = p,
    width = plot_width_in,
    height = plot_height_in,
    dpi = plot_dpi
  )
}

render_step3_dataset_plots_from_csv <- function(prefix, export_dir) {
  render_step3_dataset_gate_plot_from_csv(
    summary_csv = file.path(export_dir, paste0(prefix, '_step3_summary.csv')),
    out_png = file.path(export_dir, paste0(prefix, '_step3_gate.png')),
    title_txt = paste0('[', prefix, '] Step3 pre-fit gate summary')
  )
}

render_step3_overview_plots_from_csv <- function(export_dir, dataset_order) {
  render_step3_overview_gate_plot_from_csv(
    overview_csv = file.path(export_dir, 'step3_dataset_overview.csv'),
    out_png = file.path(export_dir, 'step3_overview_gate.png'),
    title_txt = 'Step3 pre-fit gate overview',
    dataset_order = dataset_order
  )
}

# 🔴 Execute: full Step3 pre-fit gate workflow ===============================

## 🟠 Read: Step1 summary, resolve datasets, and enforce file budget ===============================
step1_summary_df <- load_step1_summary(step1_dir)
datasets_resolved <- resolve_datasets_to_run(
  step1_summary_df = step1_summary_df,
  datasets_to_run = datasets_to_run
)

target_files <- make_target_files(datasets_resolved, export_dir)

if (length(target_files) > max_generated_files) {
  stop(
    'Planned output file count (', length(target_files),
    ') exceeds hard cap (', max_generated_files, ').'
  )
}

## 🟠 Remove: stale Step3 outputs before the current rerun ===============================
unlink(target_files[file.exists(target_files)], force = TRUE)

## 🟠 Run: dataset-wise Step3A, Step3B, overall gate, and per-dataset exports ===============================
run_step3_one_dataset <- function(prefix,
                                  step1_dir,
                                  step1_summary_df,
                                  export_dir,
                                  study_design = 'observational',
                                  arm_var = 'arm',
                                  alpha = 0.05,
                                  reps_nonzero_cure = 5000L,
                                  reps_xie_boot = 1000L,
                                  min_valid_xie_boot_prop = 0.80,
                                  seed_base = 20260306L,
                                  b_nonzerocure = NULL,
                                  receus_pi_threshold = 0.025,
                                  receus_r_threshold = 0.05,
                                  receus_candidate_dists = c('exp', 'weibull', 'gamma', 'llogis'),
                                  receus_link = 'logistic',
                                  rct_gate_rule = 'at_least_one_arm') {
  message('Running full Step3 gate for: ', prefix)
  
  bundle <- load_step1_bundle(prefix = prefix, step1_dir = step1_dir)
  check_step1_primary_endpoint(bundle = bundle, prefix = prefix)
  step1_row <- get_step1_summary_row(step1_summary_df = step1_summary_df, prefix = prefix)
  primary_df <- extract_primary_analysis_df(
    step1_bundle = bundle,
    arm_var = if (identical(study_design, 'rct')) arm_var else NULL
  )
  
  seed_nonzero_cure <- seed_base + dataset_seed_offset(prefix)
  seed_xie <- seed_base + 100000L + dataset_seed_offset(prefix)
  
  step3a_res <- run_nonzerocure_once(
    km_fit = bundle$fit$km_transition,
    reps = reps_nonzero_cure,
    seed = seed_nonzero_cure,
    alpha = alpha,
    b = b_nonzerocure
  )
  
  if (identical(study_design, 'observational')) {
    step3b_res <- run_step3b_observational(
      df = primary_df,
      km_fit = bundle$fit$km_transition,
      alpha = alpha,
      reps_xie_boot = reps_xie_boot,
      seed_xie = seed_xie,
      min_valid_prop = min_valid_xie_boot_prop
    )
  } else if (identical(study_design, 'rct')) {
    step3b_res <- run_step3b_rct(
      df = primary_df,
      arm_var = 'arm_value',
      pi_threshold = receus_pi_threshold,
      r_threshold = receus_r_threshold,
      candidate_dists = receus_candidate_dists,
      link_name = receus_link,
      gate_rule = rct_gate_rule
    )
  } else {
    stop('study_design must be either "observational" or "rct".')
  }
  
  summary_tbl <- build_step3_summary(
    step1_row = step1_row,
    study_design = study_design,
    step3a_tbl = step3a_res$summary,
    step3b_tbl = step3b_res$summary,
    alpha = alpha,
    reps_nonzero_cure = reps_nonzero_cure,
    reps_xie_boot = reps_xie_boot,
    seed_nonzero_cure = seed_nonzero_cure,
    seed_xie = seed_xie
  )
  
  write_csv_safe(
    summary_tbl,
    file.path(export_dir, paste0(prefix, '_step3_summary.csv'))
  )
  
  step3_rds <- list(
    meta = list(
      created_at = Sys.time(),
      analysis_step = 'Step3_full_pre_fit_gate',
      study_design = study_design,
      workflow_note = 'Step3 includes Step3A non-zero cure screening and Step3B sufficient follow-up / cure-appropriateness gate.',
      primary_endpoint_definition = list(
        event = 'transition (status_num == 1)',
        censoring = 'right_censoring (0) + remission (2)'
      ),
      step3a = list(
        function_name = 'hdcuremodels::nonzerocure_test',
        alpha = alpha,
        reps = reps_nonzero_cure,
        seed = seed_nonzero_cure,
        b = b_nonzerocure
      ),
      step3b = if (identical(study_design, 'observational')) {
        list(
          primary_method = 'Xie_2024_bootstrap',
          sensitivity_method = 'hdcuremodels::sufficient_fu_test',
          alpha = alpha,
          reps_xie_boot = reps_xie_boot,
          seed_xie = seed_xie
        )
      } else {
        list(
          primary_method = 'RECeUS_AIC_armwise',
          arm_var = arm_var,
          pi_threshold = receus_pi_threshold,
          r_threshold = receus_r_threshold,
          candidate_dists = receus_candidate_dists,
          link = receus_link,
          gate_rule = rct_gate_rule
        )
      },
      future_use = list(
        contains_survfit_object = TRUE,
        contains_individual_level_step1_input = TRUE,
        note = 'This RDS carries forward Step1 input and fitted KM objects so that later survival summaries can be recalculated from this file alone.'
      ),
      note = 'PNG plots are rendered from exported CSV files only.'
    ),
    input_files = list(
      step1_summary_csv = step1_summary_path(step1_dir),
      step1_bundle_rds = step1_bundle_path(prefix, step1_dir)
    ),
    step1_carryforward = build_step3_carryforward(bundle),
    screen_nonzero_cure = step3a_res$raw,
    step3b_raw = step3b_res$raw,
    summary = summary_tbl
  )
  
  saveRDS(
    step3_rds,
    file.path(export_dir, paste0(prefix, '_step3_bundle.rds'))
  )
  
  render_step3_dataset_plots_from_csv(
    prefix = prefix,
    export_dir = export_dir
  )
  
  message('Done full Step3 gate: ', prefix)
  invisible(step3_rds)
}

step3_results <- lapply(
  datasets_resolved,
  function(prefix) {
    run_step3_one_dataset(
      prefix = prefix,
      step1_dir = step1_dir,
      step1_summary_df = step1_summary_df,
      export_dir = export_dir,
      study_design = study_design,
      arm_var = arm_var,
      alpha = alpha_primary,
      reps_nonzero_cure = reps_nonzero_cure,
      reps_xie_boot = reps_xie_boot,
      min_valid_xie_boot_prop = min_valid_xie_boot_prop,
      seed_base = seed_base,
      b_nonzerocure = b_nonzerocure,
      receus_pi_threshold = receus_pi_threshold,
      receus_r_threshold = receus_r_threshold,
      receus_candidate_dists = receus_candidate_dists,
      receus_link = receus_link,
      rct_gate_rule = rct_gate_rule
    )
  }
)
names(step3_results) <- datasets_resolved

## 🟠 Collect: named screening objects and combined overview table ===============================
screen_nonzero_cure <- lapply(step3_results, function(x) x$screen_nonzero_cure)

if (identical(study_design, 'observational')) {
  screen_fu_xie <- lapply(step3_results, function(x) x$step3b_raw$xie)
  screen_fu_mz <- lapply(step3_results, function(x) x$step3b_raw$mz)
  screen_receus <- NULL
} else {
  screen_fu_xie <- NULL
  screen_fu_mz <- NULL
  screen_receus <- lapply(step3_results, function(x) x$step3b_raw)
}

overview_tbl <- dplyr::bind_rows(
  lapply(step3_results, function(x) x$summary)
) |>
  dplyr::mutate(
    file_prefix = factor(file_prefix, levels = datasets_resolved)
  ) |>
  dplyr::arrange(file_prefix) |>
  dplyr::mutate(
    file_prefix = as.character(file_prefix)
  )

write_csv_safe(
  overview_tbl,
  file.path(export_dir, 'step3_dataset_overview.csv')
)

render_step3_overview_plots_from_csv(
  export_dir = export_dir,
  dataset_order = datasets_resolved
)

## 🟠 Save: combined RDS object, session info, and manifest ===============================
step3_combined <- list(
  meta = list(
    created_at = Sys.time(),
    analysis_step = 'Step3_full_pre_fit_gate',
    study_design = study_design,
    workflow_note = 'Pass requires Step3A support and Step3B support. Borderline indicates sensitivity-only use of cure models.'
  ),
  dataset_order = datasets_resolved,
  screen_nonzero_cure = screen_nonzero_cure,
  screen_fu_xie = screen_fu_xie,
  screen_fu_mz = screen_fu_mz,
  screen_receus = screen_receus,
  results = step3_results,
  overview = overview_tbl
)

saveRDS(
  step3_combined,
  file.path(export_dir, 'step3_pre_fit_gate_all.rds')
)

capture.output(
  sessionInfo(),
  file = file.path(export_dir, 'step3_sessionInfo.txt')
)

manifest_files <- unique(c(
  list.files(export_dir, full.names = FALSE),
  'step3_manifest.csv'
))

manifest <- tibble::tibble(
  file = sort(unique(manifest_files))
)

write_csv_safe(
  manifest,
  file.path(export_dir, 'step3_manifest.csv')
)

message('Full Step3 pre-fit gate finished successfully.')
message('Generated files: ', length(target_files), ' / ', max_generated_files)