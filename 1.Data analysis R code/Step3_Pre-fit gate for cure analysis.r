# 🔴 Configure: paths and thresholds ===============================

## 🟠 Declare: editable file locations ===============================
DATA_PATH <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
EXPORT_PATH <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦3.Step3_non-zero cure fraction screening/attachments'

## 🟠 Set: reproducibility and scale options ===============================
GLOBAL_SEED <- 20260301L
DAYS_PER_YEAR <- 365.25
TIME_SCALE <- "years"

## 🟠 Set: analysis options for Step 3 ===============================
SITE_INCLUDE <- NULL
SUBGROUP_VARS <- NULL
REMISSION_AS_CENSORING <- TRUE

## 🟠 Set: Step 3A screening controls ===============================
NONZERO_REPS <- 2000L
NONZERO_ALPHA <- 0.05

## 🟠 Set: Step 3B strict follow-up controls ===============================
XIE_ALPHA <- 0.05
XIE_N_BOOT <- 500L
XIE_BOOT_MAX_ATTEMPTS_MULTIPLIER <- 10L
XIE_EPSILON_RATE <- 7 / 8

## 🟠 Set: Step 3C practical follow-up controls ===============================
RUN_PRACTICAL_SUFFICIENCY <- TRUE
PRACTICAL_METHOD <- "sg"
PRACTICAL_EPS <- 0.01
PRACTICAL_ALPHA <- 0.05
PRACTICAL_N_BOOT <- 1000L
PRACTICAL_TAU_BUFFER_PROP <- 0.05
PRACTICAL_TAU_OVERRIDE_YEARS <- NULL
AUTO_INSTALL_GITHUB_PACKAGES <- TRUE

## 🟠 Set: export file names ===============================
FILE_ANALYSIS_UNITS_CSV <- "step3_analysis_units.csv"
FILE_NONZERO_CSV <- "step3_nonzero_cure_results.csv"
FILE_STRICT_CSV <- "step3_strict_followup_xie_results.csv"
FILE_PRACTICAL_CSV <- "step3_practical_followup_results.csv"
FILE_SUMMARY_CSV <- "step3_summary_results.csv"
FILE_FULL_RDS <- "step3_full_results.rds"

# 🔴 Attach: package dependencies ===============================

## 🟠 Define: package installation helpers ===============================
install_if_missing <- function(pkgs, repos = "https://cloud.r-project.org") {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0L) {
    install.packages(missing_pkgs, repos = repos, dependencies = TRUE)
  }
  invisible(NULL)
}

install_github_if_missing <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE) && isTRUE(AUTO_INSTALL_GITHUB_PACKAGES)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org", dependencies = TRUE)
    }
    remotes::install_github(repo, dependencies = TRUE, upgrade = "never")
  }
  invisible(requireNamespace(pkg, quietly = TRUE))
}

## 🟠 Load: required libraries ===============================
install_if_missing(c("survival", "hdcuremodels", "readr", "dplyr", "purrr", "tibble", "stringr"))

suppressPackageStartupMessages({
  library(survival)
  library(hdcuremodels)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(stringr)
})

set.seed(GLOBAL_SEED)
options(stringsAsFactors = FALSE)

# 🔴 Prepare: export directory and runtime log ===============================

## 🟠 Create: target export folder ===============================
if (!dir.exists(EXPORT_PATH)) {
  dir.create(EXPORT_PATH, recursive = TRUE, showWarnings = FALSE)
}

## 🟠 Initialize: runtime containers ===============================
runtime_log <- list(
  practical_pkg_available = FALSE,
  practical_pkg_name = "cureSFUTest",
  practical_pkg_repo = "tp-yuen/cureSFUTest",
  notes = character(0)
)

# 🔴 Define: helper utilities ===============================

## 🟠 Create: generic helpers ===============================
`%||%` <- function(x, y) if (is.null(x)) y else x

safe_scalar <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  as.numeric(x[[1]])
}

safe_character <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_character_)
  as.character(x[[1]])
}

append_note <- function(x, note) {
  if (is.null(x) || length(x) == 0L || is.na(x) || note == "") {
    return(note)
  }
  paste(x, note, sep = " | ")
}

assert_required_columns <- function(dat, required_cols) {
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop("필수 컬럼이 없습니다: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

resolve_status_num <- function(dat) {
  if ("status_num" %in% names(dat)) {
    out <- suppressWarnings(as.integer(dat$status_num))
    return(out)
  }
  
  if ("status" %in% names(dat)) {
    status_chr <- tolower(trimws(as.character(dat$status)))
    out <- dplyr::case_when(
      status_chr %in% c("right_censoring", "right censoring", "censor", "censored") ~ 0L,
      status_chr %in% c("transition", "event", "failure") ~ 1L,
      status_chr %in% c("remission") ~ 2L,
      TRUE ~ NA_integer_
    )
    return(out)
  }
  
  stop("status_num 또는 status 컬럼이 필요합니다.")
}

resolve_time_years <- function(dat) {
  if (!"days_followup" %in% names(dat)) {
    stop("days_followup 컬럼이 필요합니다.")
  }
  as.numeric(dat$days_followup) / DAYS_PER_YEAR
}

sanitize_level <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9_\\-]+", "_")
  x <- stringr::str_replace_all(x, "_{2,}", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  ifelse(nchar(x) == 0L, "NA", x)
}

# 🔴 Define: data engineering helpers ===============================

## 🟠 Transform: raw data into transition-only endpoint ===============================
build_analysis_data <- function(dat) {
  dat2 <- dat %>%
    mutate(
      site = as.character(site),
      id = as.character(id),
      subject_uid = paste(site, id, sep = "__"),
      status_num_resolved = resolve_status_num(cur_data_all()),
      time_years = resolve_time_years(cur_data_all()),
      event_transition = ifelse(status_num_resolved == 1L, 1L, 0L),
      censor_reason = case_when(
        status_num_resolved == 1L ~ "transition",
        status_num_resolved == 2L ~ "remission_as_censor",
        status_num_resolved == 0L ~ "right_censoring",
        TRUE ~ "unknown"
      )
    ) %>%
    filter(!is.na(site), !is.na(id), !is.na(time_years), !is.na(event_transition)) %>%
    filter(time_years >= 0)
  
  if (!is.null(SITE_INCLUDE)) {
    dat2 <- dat2 %>% filter(site %in% SITE_INCLUDE)
  }
  
  if (nrow(dat2) == 0L) {
    stop("전처리 후 남은 데이터가 없습니다.")
  }
  
  dat2
}

## 🟠 Construct: analysis units for merged and site-specific branches ===============================
build_analysis_units <- function(dat, subgroup_vars = NULL) {
  units <- list()
  
  units[[length(units) + 1L]] <- list(
    analysis_id = "MERGED__overall",
    dataset_label = "MERGED",
    subgroup_var = "overall",
    subgroup_level = "overall",
    data = dat
  )
  
  site_levels <- sort(unique(dat$site))
  for (s in site_levels) {
    units[[length(units) + 1L]] <- list(
      analysis_id = paste0(sanitize_level(s), "__overall"),
      dataset_label = as.character(s),
      subgroup_var = "overall",
      subgroup_level = "overall",
      data = dat %>% filter(site == s)
    )
  }
  
  if (!is.null(subgroup_vars) && length(subgroup_vars) > 0L) {
    for (sv in subgroup_vars) {
      if (!sv %in% names(dat)) next
      
      subgroup_levels <- unique(dat[[sv]])
      subgroup_levels <- subgroup_levels[!is.na(subgroup_levels)]
      
      for (lvl in subgroup_levels) {
        units[[length(units) + 1L]] <- list(
          analysis_id = paste0("MERGED__", sanitize_level(sv), "__", sanitize_level(lvl)),
          dataset_label = "MERGED",
          subgroup_var = sv,
          subgroup_level = as.character(lvl),
          data = dat %>% filter(.data[[sv]] == lvl)
        )
        
        for (s in site_levels) {
          tmp <- dat %>% filter(site == s, .data[[sv]] == lvl)
          if (nrow(tmp) == 0L) next
          units[[length(units) + 1L]] <- list(
            analysis_id = paste0(sanitize_level(s), "__", sanitize_level(sv), "__", sanitize_level(lvl)),
            dataset_label = as.character(s),
            subgroup_var = sv,
            subgroup_level = as.character(lvl),
            data = tmp
          )
        }
      }
    }
  }
  
  units
}

## 🟠 Summarize: analysis unit basics ===============================
summarise_unit <- function(unit) {
  d <- unit$data
  tibble(
    analysis_id = unit$analysis_id,
    dataset_label = unit$dataset_label,
    subgroup_var = unit$subgroup_var,
    subgroup_level = unit$subgroup_level,
    n = nrow(d),
    n_event = sum(d$event_transition == 1L, na.rm = TRUE),
    n_censor = sum(d$event_transition == 0L, na.rm = TRUE),
    max_time_years = suppressWarnings(max(d$time_years, na.rm = TRUE)),
    last_event_time_years = ifelse(sum(d$event_transition == 1L, na.rm = TRUE) > 0L,
                                   max(d$time_years[d$event_transition == 1L], na.rm = TRUE),
                                   NA_real_),
    last_censor_time_years = ifelse(sum(d$event_transition == 0L, na.rm = TRUE) > 0L,
                                    max(d$time_years[d$event_transition == 0L], na.rm = TRUE),
                                    NA_real_)
  )
}

# 🔴 Define: Kaplan-Meier helpers ===============================

## 🟠 Fit: plain Kaplan-Meier for one analysis unit ===============================
fit_km_unit <- function(unit) {
  d <- unit$data
  if (nrow(d) == 0L) return(NULL)
  survfit(Surv(time_years, event_transition) ~ 1, data = d)
}

# 🔴 Define: Step 3A helpers ===============================

## 🟠 Run: non-zero cure screening on survfit object ===============================
run_nonzero_cure_screen <- function(km_fit, alpha = NONZERO_ALPHA, reps = NONZERO_REPS, seed = GLOBAL_SEED) {
  if (is.null(km_fit)) {
    return(list(
      ok = FALSE,
      result = NULL,
      note = "km_fit_is_null"
    ))
  }
  
  out <- tryCatch(
    hdcuremodels::nonzerocure_test(
      object = km_fit,
      reps = reps,
      seed = seed,
      plot = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(out, "error")) {
    return(list(
      ok = FALSE,
      result = NULL,
      note = paste0("nonzero_test_error: ", conditionMessage(out))
    ))
  }
  
  result_tbl <- tibble(
    proportion_susceptible = safe_scalar(out$proportion_susceptible),
    proportion_cured = safe_scalar(out$proportion_cured),
    p_value = safe_scalar(out$p_value),
    time_95_percent_of_events = safe_scalar(out$time_95_percent_of_events),
    alpha = alpha,
    reject_h0_no_cure = ifelse(!is.na(safe_scalar(out$p_value)) && safe_scalar(out$p_value) < alpha, TRUE, FALSE),
    interpretation = case_when(
      is.na(safe_scalar(out$p_value)) ~ "not_run",
      safe_scalar(out$p_value) < alpha ~ "support_nonzero_cure_signal",
      TRUE ~ "insufficient_evidence_for_nonzero_cure_signal"
    )
  )
  
  list(ok = TRUE, result = out, table = result_tbl, note = NA_character_)
}

# 🔴 Define: Step 3B helpers for strict follow-up ===============================

## 🟠 Implement: Xie-type Kaplan-Meier cdf estimator ===============================
xie_Fhat <- function(y, delta, t) {
  if (length(y) == 0L || all(is.na(y)) || is.na(t)) return(NA_real_)
  y <- as.numeric(y)
  delta <- as.integer(delta)
  
  if (t < min(y, na.rm = TRUE)) return(0)
  
  ord <- order(y)
  y <- y[ord]
  delta <- delta[ord]
  n <- length(y)
  
  surv_mult <- 1 - delta / (n:1)
  Fhat_vals <- 1 - cumprod(surv_mult)
  idx <- max(which(y <= t))
  Fhat_vals[idx]
}

## 🟠 Choose: epsilon for strict follow-up test ===============================
choose_xie_epsilon <- function(y, delta, phi = XIE_EPSILON_RATE) {
  y <- as.numeric(y)
  delta <- as.integer(delta)
  
  if (sum(delta == 1L, na.rm = TRUE) == 0L) return(NA_real_)
  
  t_max <- max(y, na.rm = TRUE)
  tf_max <- max(y[delta == 1L], na.rm = TRUE)
  gap <- t_max - tf_max
  
  if (is.na(t_max) || is.na(tf_max)) return(NA_real_)
  
  if (2 * gap >= t_max) {
    t_max
  } else {
    2 * gap + phi * (t_max - 2 * gap)
  }
}

## 🟠 Compute: strict follow-up statistic ===============================
xie_statistic <- function(y, delta, epsilon) {
  y <- as.numeric(y)
  delta <- as.integer(delta)
  ep <- as.numeric(epsilon)
  
  if (anyNA(c(ep)) || length(y) == 0L || sum(delta == 1L, na.rm = TRUE) == 0L) {
    return(NA_real_)
  }
  
  t_max <- max(y, na.rm = TRUE)
  tf_max <- max(y[delta == 1L], na.rm = TRUE)
  
  if (ep <= (t_max - tf_max) || ep > t_max) {
    return(NA_real_)
  }
  
  F0 <- xie_Fhat(y, delta, t_max)
  
  pGhat <- function(t) {
    num <- (xie_Fhat(y, delta, t_max - t / 2) - xie_Fhat(y, delta, t_max - t))^2
    den <- 2 * xie_Fhat(y, delta, t_max - t / 2) - xie_Fhat(y, delta, t_max - t) - xie_Fhat(y, delta, t_max)
    xie_Fhat(y, delta, t_max - t) + num / den
  }
  
  pGhat_value <- suppressWarnings(pGhat(ep))
  
  if (is.na(pGhat_value)) return(0)
  if (pGhat_value > 1) return(1 - F0)
  if (pGhat_value < F0) return(0)
  
  pGhat_value - F0
}

## 🟠 Bootstrap: strict follow-up critical values ===============================
xie_bootstrap_centered <- function(y, delta, epsilon, n_boot = XIE_N_BOOT,
                                   max_attempts_multiplier = XIE_BOOT_MAX_ATTEMPTS_MULTIPLIER,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  y <- as.numeric(y)
  delta <- as.integer(delta)
  n <- length(y)
  t_obs <- xie_statistic(y, delta, epsilon)
  
  if (is.na(t_obs)) {
    return(list(centered = numeric(0), statistic = NA_real_, n_success = 0L, n_attempt = 0L))
  }
  
  target <- as.integer(n_boot)
  max_attempts <- as.integer(max(1000L, target * max_attempts_multiplier))
  centered <- numeric(0)
  attempts <- 0L
  
  while (length(centered) < target && attempts < max_attempts) {
    attempts <- attempts + 1L
    idx <- sample.int(n = n, size = n, replace = TRUE)
    y_b <- y[idx]
    delta_b <- delta[idx]
    
    if (delta_b[which.max(y_b)] == 1L) next
    
    t_boot <- xie_statistic(y_b, delta_b, epsilon)
    if (is.na(t_boot) || is.null(t_boot)) next
    
    centered <- c(centered, t_boot - t_obs)
  }
  
  list(
    centered = centered,
    statistic = t_obs,
    n_success = length(centered),
    n_attempt = attempts
  )
}

## 🟠 Run: strict follow-up test for one analysis unit ===============================
run_xie_strict_followup <- function(df_unit,
                                    alpha = XIE_ALPHA,
                                    n_boot = XIE_N_BOOT,
                                    phi = XIE_EPSILON_RATE,
                                    seed = GLOBAL_SEED) {
  y <- df_unit$time_years
  delta <- df_unit$event_transition
  
  if (length(y) == 0L) {
    return(list(ok = FALSE, result = NULL, note = "empty_unit"))
  }
  
  if (sum(delta == 1L, na.rm = TRUE) == 0L) {
    return(list(ok = FALSE, result = NULL, note = "no_event_observed"))
  }
  
  epsilon <- choose_xie_epsilon(y, delta, phi = phi)
  boot_obj <- xie_bootstrap_centered(y, delta, epsilon, n_boot = n_boot, seed = seed)
  
  if (length(boot_obj$centered) < 20L || is.na(boot_obj$statistic)) {
    return(list(
      ok = FALSE,
      result = NULL,
      note = paste0("insufficient_valid_bootstrap_samples: ", boot_obj$n_success)
    ))
  }
  
  crit_val <- as.numeric(stats::quantile(boot_obj$centered, probs = 1 - alpha, na.rm = TRUE, type = 7))
  p_value <- mean(boot_obj$centered >= boot_obj$statistic, na.rm = TRUE)
  
  result_tbl <- tibble(
    epsilon = epsilon,
    statistic = boot_obj$statistic,
    crit_val = crit_val,
    p_value = p_value,
    alpha = alpha,
    n_boot_target = n_boot,
    n_boot_success = boot_obj$n_success,
    n_boot_attempts = boot_obj$n_attempt,
    reject_h0_sufficient_followup = ifelse(!is.na(p_value) && p_value < alpha, TRUE, FALSE),
    interpretation = case_when(
      is.na(p_value) ~ "not_run",
      p_value < alpha ~ "evidence_of_insufficient_followup",
      TRUE ~ "no_evidence_against_sufficient_followup"
    )
  )
  
  list(
    ok = TRUE,
    result = list(
      epsilon = epsilon,
      statistic = boot_obj$statistic,
      centered_bootstrap = boot_obj$centered,
      crit_val = crit_val,
      p_value = p_value,
      alpha = alpha,
      n_boot_target = n_boot,
      n_boot_success = boot_obj$n_success,
      n_boot_attempts = boot_obj$n_attempt
    ),
    table = result_tbl,
    note = NA_character_
  )
}

# 🔴 Define: Step 3C helpers for practical follow-up ===============================

## 🟠 Resolve: practical sufficiency package availability ===============================
ensure_practical_package <- function() {
  if (requireNamespace("cureSFUTest", quietly = TRUE)) {
    runtime_log$practical_pkg_available <<- TRUE
    return(TRUE)
  }
  
  ok <- FALSE
  if (isTRUE(AUTO_INSTALL_GITHUB_PACKAGES)) {
    ok <- tryCatch(
      {
        install_github_if_missing("cureSFUTest", "tp-yuen/cureSFUTest")
        requireNamespace("cureSFUTest", quietly = TRUE)
      },
      error = function(e) FALSE
    )
  }
  
  runtime_log$practical_pkg_available <<- isTRUE(ok)
  isTRUE(ok)
}

## 🟠 Resolve: practical tau value ===============================
resolve_practical_tau <- function(df_unit, analysis_id, dataset_label) {
  y_max <- max(df_unit$time_years, na.rm = TRUE)
  
  if (is.null(PRACTICAL_TAU_OVERRIDE_YEARS)) {
    tau_g <- y_max
    tau <- if (tau_g > 0) tau_g * (1 + PRACTICAL_TAU_BUFFER_PROP) else 1e-06
    return(list(tau = tau, tau_g = tau_g))
  }
  
  if (length(PRACTICAL_TAU_OVERRIDE_YEARS) == 1L && is.numeric(PRACTICAL_TAU_OVERRIDE_YEARS)) {
    tau_g <- y_max
    tau <- as.numeric(PRACTICAL_TAU_OVERRIDE_YEARS)
    return(list(tau = tau, tau_g = tau_g))
  }
  
  if (!is.null(names(PRACTICAL_TAU_OVERRIDE_YEARS))) {
    if (analysis_id %in% names(PRACTICAL_TAU_OVERRIDE_YEARS)) {
      tau_g <- y_max
      tau <- as.numeric(PRACTICAL_TAU_OVERRIDE_YEARS[[analysis_id]])
      return(list(tau = tau, tau_g = tau_g))
    }
    if (dataset_label %in% names(PRACTICAL_TAU_OVERRIDE_YEARS)) {
      tau_g <- y_max
      tau <- as.numeric(PRACTICAL_TAU_OVERRIDE_YEARS[[dataset_label]])
      return(list(tau = tau, tau_g = tau_g))
    }
  }
  
  tau_g <- y_max
  tau <- if (tau_g > 0) tau_g * (1 + PRACTICAL_TAU_BUFFER_PROP) else 1e-06
  list(tau = tau, tau_g = tau_g)
}

## 🟠 Run: practical sufficient follow-up test ===============================
run_practical_followup <- function(df_unit,
                                   analysis_id,
                                   dataset_label,
                                   alpha = PRACTICAL_ALPHA,
                                   eps = PRACTICAL_EPS,
                                   n_boot = PRACTICAL_N_BOOT,
                                   method = PRACTICAL_METHOD) {
  if (!isTRUE(RUN_PRACTICAL_SUFFICIENCY)) {
    return(list(ok = FALSE, result = NULL, note = "practical_test_disabled"))
  }
  
  if (!ensure_practical_package()) {
    return(list(ok = FALSE, result = NULL, note = "cureSFUTest_unavailable"))
  }
  
  y <- as.numeric(df_unit$time_years)
  delta <- as.integer(df_unit$event_transition)
  
  if (length(y) == 0L) {
    return(list(ok = FALSE, result = NULL, note = "empty_unit"))
  }
  
  tau_obj <- resolve_practical_tau(df_unit, analysis_id, dataset_label)
  
  out <- tryCatch(
    cureSFUTest::sfu.test(
      y = y,
      delta = delta,
      tau = tau_obj$tau,
      tau.g = tau_obj$tau_g,
      eps = eps,
      method = method,
      alpha = alpha,
      n.boot = n_boot
    ),
    error = function(e) e
  )
  
  if (inherits(out, "error")) {
    return(list(ok = FALSE, result = NULL, note = paste0("practical_test_error: ", conditionMessage(out))))
  }
  
  result_tbl <- tibble(
    statistic = safe_scalar(out$statistic),
    p_value = safe_scalar(out$p.value),
    crit_val = safe_scalar(out$crit.val),
    method = safe_character(out$method),
    alpha = safe_scalar(out$alpha),
    tau_g = safe_scalar(out$tau.g),
    tau = safe_scalar(out$tau),
    eps = safe_scalar(out$eps),
    reject_h0_practical_insufficient_followup = ifelse(!is.na(safe_scalar(out$p.value)) && safe_scalar(out$p.value) < alpha, TRUE, FALSE),
    interpretation = case_when(
      is.na(safe_scalar(out$p.value)) ~ "not_run",
      safe_scalar(out$p.value) < alpha ~ "support_practical_sufficient_followup",
      TRUE ~ "insufficient_evidence_for_practical_sufficient_followup"
    )
  )
  
  list(ok = TRUE, result = out, table = result_tbl, note = NA_character_)
}

# 🔴 Read: raw input and derive transition endpoint ===============================

## 🟠 Import: merged cohort dataset ===============================
raw_dat <- readr::read_csv(DATA_PATH, show_col_types = FALSE)

## 🟠 Validate: required columns before transformation ===============================
assert_required_columns(raw_dat, c("id", "site", "days_followup"))

## 🟠 Build: transition-only analysis data ===============================
analysis_dat <- build_analysis_data(raw_dat)

# 🔴 Construct: analysis units and Kaplan-Meier objects ===============================

## 🟠 Split: merged and site-specific branches ===============================
analysis_units <- build_analysis_units(analysis_dat, subgroup_vars = SUBGROUP_VARS)

## 🟠 Fit: plain Kaplan-Meier objects for all units ===============================
unit_results <- purrr::map(
  analysis_units,
  function(unit) {
    km_fit <- fit_km_unit(unit)
    unit_info <- summarise_unit(unit)
    list(
      meta = unit,
      info = unit_info,
      km_fit = km_fit
    )
  }
)

## 🟠 Assemble: analysis unit registry table ===============================
analysis_units_tbl <- purrr::map_dfr(unit_results, "info")

# 🔴 Run: Step 3A non-zero cure screening ===============================

## 🟠 Execute: screening across analysis units ===============================
step3a_list <- purrr::imap(
  unit_results,
  function(x, idx) {
    meta <- x$meta
    km_fit <- x$km_fit
    info <- x$info
    
    res <- run_nonzero_cure_screen(
      km_fit = km_fit,
      alpha = NONZERO_ALPHA,
      reps = NONZERO_REPS,
      seed = GLOBAL_SEED + idx
    )
    
    tbl <- info %>%
      mutate(
        step = "3A",
        test_name = "nonzero_cure_screen_hdcuremodels",
        test_ok = res$ok,
        note = res$note
      )
    
    if (isTRUE(res$ok)) {
      tbl <- bind_cols(tbl, res$table)
    } else {
      tbl <- tbl %>%
        mutate(
          proportion_susceptible = NA_real_,
          proportion_cured = NA_real_,
          p_value = NA_real_,
          time_95_percent_of_events = NA_real_,
          alpha = NONZERO_ALPHA,
          reject_h0_no_cure = NA,
          interpretation = "not_run"
        )
    }
    
    list(meta = meta, raw = res$result, table = tbl)
  }
)

step3a_tbl <- purrr::map_dfr(step3a_list, "table")

# 🔴 Run: Step 3B strict follow-up via Xie-type implementation ===============================

## 🟠 Execute: strict follow-up tests across analysis units ===============================
step3b_list <- purrr::imap(
  unit_results,
  function(x, idx) {
    meta <- x$meta
    d <- meta$data
    info <- x$info
    
    res <- run_xie_strict_followup(
      df_unit = d,
      alpha = XIE_ALPHA,
      n_boot = XIE_N_BOOT,
      phi = XIE_EPSILON_RATE,
      seed = GLOBAL_SEED + 1000L + idx
    )
    
    tbl <- info %>%
      mutate(
        step = "3B",
        test_name = "strict_followup_xie_extremes",
        test_ok = res$ok,
        note = res$note
      )
    
    if (isTRUE(res$ok)) {
      tbl <- bind_cols(tbl, res$table)
    } else {
      tbl <- tbl %>%
        mutate(
          epsilon = NA_real_,
          statistic = NA_real_,
          crit_val = NA_real_,
          p_value = NA_real_,
          alpha = XIE_ALPHA,
          n_boot_target = XIE_N_BOOT,
          n_boot_success = NA_integer_,
          n_boot_attempts = NA_integer_,
          reject_h0_sufficient_followup = NA,
          interpretation = "not_run"
        )
    }
    
    list(meta = meta, raw = res$result, table = tbl)
  }
)

step3b_tbl <- purrr::map_dfr(step3b_list, "table")

# 🔴 Run: Step 3C practical follow-up via Yuen-Musta ===============================

## 🟠 Execute: practical follow-up tests across analysis units ===============================
step3c_list <- purrr::imap(
  unit_results,
  function(x, idx) {
    meta <- x$meta
    d <- meta$data
    info <- x$info
    
    res <- run_practical_followup(
      df_unit = d,
      analysis_id = meta$analysis_id,
      dataset_label = meta$dataset_label,
      alpha = PRACTICAL_ALPHA,
      eps = PRACTICAL_EPS,
      n_boot = PRACTICAL_N_BOOT,
      method = PRACTICAL_METHOD
    )
    
    tbl <- info %>%
      mutate(
        step = "3C",
        test_name = "practical_followup_yuen_musta",
        test_ok = res$ok,
        note = res$note
      )
    
    if (isTRUE(res$ok)) {
      tbl <- bind_cols(tbl, res$table)
    } else {
      tbl <- tbl %>%
        mutate(
          statistic = NA_real_,
          p_value = NA_real_,
          crit_val = NA_real_,
          method = PRACTICAL_METHOD,
          alpha = PRACTICAL_ALPHA,
          tau_g = NA_real_,
          tau = NA_real_,
          eps = PRACTICAL_EPS,
          reject_h0_practical_insufficient_followup = NA,
          interpretation = "not_run"
        )
    }
    
    list(meta = meta, raw = res$result, table = tbl)
  }
)

step3c_tbl <- purrr::map_dfr(step3c_list, "table")

# 🔴 Summarize: Step 3 outputs into a single table ===============================

## 🟠 Join: Step 3A, 3B, and 3C summaries ===============================
step3_summary_tbl <- analysis_units_tbl %>%
  select(analysis_id, dataset_label, subgroup_var, subgroup_level, n, n_event, n_censor,
         max_time_years, last_event_time_years, last_censor_time_years) %>%
  left_join(
    step3a_tbl %>%
      select(
        analysis_id,
        nz_proportion_susceptible = proportion_susceptible,
        nz_proportion_cured = proportion_cured,
        nz_p_value = p_value,
        nz_time_95_percent_of_events = time_95_percent_of_events,
        nz_reject_h0_no_cure = reject_h0_no_cure,
        nz_interpretation = interpretation,
        nz_note = note
      ),
    by = "analysis_id"
  ) %>%
  left_join(
    step3b_tbl %>%
      select(
        analysis_id,
        strict_epsilon = epsilon,
        strict_statistic = statistic,
        strict_crit_val = crit_val,
        strict_p_value = p_value,
        strict_reject_h0_sufficient_followup = reject_h0_sufficient_followup,
        strict_interpretation = interpretation,
        strict_note = note
      ),
    by = "analysis_id"
  ) %>%
  left_join(
    step3c_tbl %>%
      select(
        analysis_id,
        practical_statistic = statistic,
        practical_p_value = p_value,
        practical_crit_val = crit_val,
        practical_method = method,
        practical_tau_g = tau_g,
        practical_tau = tau,
        practical_eps = eps,
        practical_reject_h0_insufficient_followup = reject_h0_practical_insufficient_followup,
        practical_interpretation = interpretation,
        practical_note = note
      ),
    by = "analysis_id"
  ) %>%
  mutate(
    step3_nonzero_signal = case_when(
      is.na(nz_reject_h0_no_cure) ~ "not_run",
      nz_reject_h0_no_cure ~ "present",
      TRUE ~ "not_supported"
    ),
    step3_strict_followup = case_when(
      is.na(strict_reject_h0_sufficient_followup) ~ "not_run",
      strict_reject_h0_sufficient_followup ~ "insufficient_followup_signal",
      TRUE ~ "no_evidence_against_sufficient_followup"
    ),
    step3_practical_followup = case_when(
      is.na(practical_reject_h0_insufficient_followup) ~ "not_run",
      practical_reject_h0_insufficient_followup ~ "support_practical_sufficient_followup",
      TRUE ~ "insufficient_evidence_for_practical_sufficient_followup"
    ),
    step3_profile = case_when(
      step3_nonzero_signal == "present" &
        step3_strict_followup == "no_evidence_against_sufficient_followup" &
        step3_practical_followup == "support_practical_sufficient_followup" ~ "supportive_on_all_three",
      step3_nonzero_signal == "present" &
        step3_strict_followup == "no_evidence_against_sufficient_followup" &
        step3_practical_followup != "support_practical_sufficient_followup" ~ "nonzero_plus_strict_only",
      step3_nonzero_signal == "present" &
        step3_strict_followup != "no_evidence_against_sufficient_followup" &
        step3_practical_followup == "support_practical_sufficient_followup" ~ "nonzero_plus_practical_only",
      step3_nonzero_signal == "present" ~ "nonzero_only_or_followup_mixed",
      TRUE ~ "weak_or_mixed"
    )
  )

# 🔴 Export: csv summaries and full rds object ===============================

## 🟠 Write: step-wise csv outputs ===============================
readr::write_csv(analysis_units_tbl, file.path(EXPORT_PATH, FILE_ANALYSIS_UNITS_CSV))
readr::write_csv(step3a_tbl, file.path(EXPORT_PATH, FILE_NONZERO_CSV))
readr::write_csv(step3b_tbl, file.path(EXPORT_PATH, FILE_STRICT_CSV))
readr::write_csv(step3c_tbl, file.path(EXPORT_PATH, FILE_PRACTICAL_CSV))
readr::write_csv(step3_summary_tbl, file.path(EXPORT_PATH, FILE_SUMMARY_CSV))

## 🟠 Save: full reusable results object to rds ===============================
step3_full_results <- list(
  config = list(
    DATA_PATH = DATA_PATH,
    EXPORT_PATH = EXPORT_PATH,
    GLOBAL_SEED = GLOBAL_SEED,
    DAYS_PER_YEAR = DAYS_PER_YEAR,
    TIME_SCALE = TIME_SCALE,
    SITE_INCLUDE = SITE_INCLUDE,
    SUBGROUP_VARS = SUBGROUP_VARS,
    NONZERO_REPS = NONZERO_REPS,
    NONZERO_ALPHA = NONZERO_ALPHA,
    XIE_ALPHA = XIE_ALPHA,
    XIE_N_BOOT = XIE_N_BOOT,
    XIE_BOOT_MAX_ATTEMPTS_MULTIPLIER = XIE_BOOT_MAX_ATTEMPTS_MULTIPLIER,
    XIE_EPSILON_RATE = XIE_EPSILON_RATE,
    RUN_PRACTICAL_SUFFICIENCY = RUN_PRACTICAL_SUFFICIENCY,
    PRACTICAL_METHOD = PRACTICAL_METHOD,
    PRACTICAL_EPS = PRACTICAL_EPS,
    PRACTICAL_ALPHA = PRACTICAL_ALPHA,
    PRACTICAL_N_BOOT = PRACTICAL_N_BOOT,
    PRACTICAL_TAU_BUFFER_PROP = PRACTICAL_TAU_BUFFER_PROP,
    PRACTICAL_TAU_OVERRIDE_YEARS = PRACTICAL_TAU_OVERRIDE_YEARS
  ),
  runtime_log = runtime_log,
  raw_data = raw_dat,
  analysis_data = analysis_dat,
  analysis_units = purrr::map(analysis_units, ~ .x[c("analysis_id", "dataset_label", "subgroup_var", "subgroup_level")]),
  km_fits = purrr::set_names(purrr::map(unit_results, "km_fit"), purrr::map_chr(unit_results, ~ .x$meta$analysis_id)),
  step3a_raw = purrr::set_names(purrr::map(step3a_list, "raw"), purrr::map_chr(step3a_list, ~ .x$meta$analysis_id)),
  step3b_raw = purrr::set_names(purrr::map(step3b_list, "raw"), purrr::map_chr(step3b_list, ~ .x$meta$analysis_id)),
  step3c_raw = purrr::set_names(purrr::map(step3c_list, "raw"), purrr::map_chr(step3c_list, ~ .x$meta$analysis_id)),
  tables = list(
    analysis_units = analysis_units_tbl,
    step3a = step3a_tbl,
    step3b = step3b_tbl,
    step3c = step3c_tbl,
    step3_summary = step3_summary_tbl
  ),
  session_info = utils::capture.output(sessionInfo())
)

saveRDS(step3_full_results, file.path(EXPORT_PATH, FILE_FULL_RDS))

# 🔴 Finish: console summary for quick verification ===============================

## 🟠 Print: exported file paths and row counts ===============================
message("Step 3 완료")
message("분석 단위 수: ", nrow(analysis_units_tbl))
message("CSV export: ", file.path(EXPORT_PATH, FILE_ANALYSIS_UNITS_CSV))
message("CSV export: ", file.path(EXPORT_PATH, FILE_NONZERO_CSV))
message("CSV export: ", file.path(EXPORT_PATH, FILE_STRICT_CSV))
message("CSV export: ", file.path(EXPORT_PATH, FILE_PRACTICAL_CSV))
message("CSV export: ", file.path(EXPORT_PATH, FILE_SUMMARY_CSV))
message("RDS export: ", file.path(EXPORT_PATH, FILE_FULL_RDS))