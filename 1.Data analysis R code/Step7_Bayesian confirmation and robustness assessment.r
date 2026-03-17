# 🔴 Configure: paths and runtime ===============================
data_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦7.Step7_베이지안 MCM 방법/🟩2.메인모델링 결과/attachments'

auto_install_packages <- FALSE
overwrite_existing_exports <- TRUE
save_stanfit_inside_master_rds <- TRUE

prediction_years <- 1:10
selected_contrast_years <- c(2, 5, 10)
calibration_years <- c(2, 5, 10)

active_bundle_ids <- c(
  "primary_family",
  "sa1_robust",
  "sa2_shape_off",
  "sa3_no_external",
  "sa4_age_support",
  "sa5_ltight",
  "sa5_lwide",
  "sa7b_remission_ic",
  "sa8_delta_narrow",
  "sa8_delta_wide",
  "sa9_anchor_alt",
  "sa10_betaMale_direct",
  "sa11_int_tight",
  "sa11_int_wide",
  "sa12_noncure_benchmark"
)

run_remission_competing_risk <- TRUE

sampling_seed <- 20250312L
chains_main <- 4L
warmup_main <- 2000L
iter_sampling_main <- 2000L
adapt_delta_main <- 0.95
max_treedepth_main <- 15L

chains_refit <- 4L
warmup_refit <- 2000L
iter_sampling_refit <- 2000L
adapt_delta_refit <- 0.99
max_treedepth_refit <- 18L

prediction_draws_max <- 1000L
ess_threshold <- 400
tail_fit_absdiff_threshold <- 0.05

if (!dir.exists(export_path)) {
  dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
}

# 🔴 Initialize: packages and options ===============================
options(stringsAsFactors = FALSE)
Sys.setenv(TZ = "UTC")

required_packages <- c(
  "rstan", "survival", "cmprsk", "loo", "posterior", "timeROC",
  "dplyr", "tibble", "tidyr", "purrr", "readr", "stringr", "data.table"
)

load_or_install_packages <- function(pkgs, auto_install = FALSE) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0L) {
    if (!auto_install) {
      stop(
        "다음 패키지가 설치되어 있지 않습니다: ",
        paste(missing_pkgs, collapse = ", "),
        ". auto_install_packages <- TRUE 로 바꾸거나 수동 설치 후 다시 실행하세요."
      )
    }
    install.packages(missing_pkgs, dependencies = TRUE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

load_or_install_packages(required_packages, auto_install = auto_install_packages)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(parallel::detectCores(), chains_main))

# 🔴 Define: shared objects and constants ===============================
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

bind_rows_or_template <- function(lst, template = NULL) {
  keep <- lst[!vapply(lst, is.null, logical(1))]
  if (length(keep) == 0L) {
    if (is.null(template)) return(tibble::tibble())
    return(template[0, , drop = FALSE])
  }
  dplyr::bind_rows(keep)
}

family_lookup <- tibble::tibble(
  family_id = 1L:4L,
  family_name = c("exponential", "weibull", "loglogistic", "lognormal")
)

bundle_latency_map <- list(
  primary_family = list(
    PNU = c("B-L0", "B-L1"),
    SNU = c("B-L0", "B-L1"),
    merged = c("B-L0S0", "B-L1S0", "B-L0S1", "B-L1S1")
  ),
  sa1_robust = list(
    PNU = c("B-L0"),
    SNU = c("B-L0"),
    merged = c("B-L0S0", "B-L0S1")
  ),
  sa2_shape_off = list(
    PNU = c("B-L0"),
    SNU = c("B-L0"),
    merged = c("B-L0S0", "B-L0S1")
  ),
  sa3_no_external = list(
    PNU = c("B-L0"),
    SNU = c("B-L0"),
    merged = c("B-L0S0", "B-L0S1")
  ),
  sa4_age_support = list(
    PNU = c("B-L0"),
    SNU = c("B-L0"),
    merged = c("B-L0S0", "B-L0S1")
  ),
  sa5_ltight = list(
    merged = c("B-L0S0", "B-L1S0", "B-L0S1", "B-L1S1")
  ),
  sa5_lwide = list(
    merged = c("B-L0S0", "B-L1S0", "B-L0S1", "B-L1S1")
  ),
  sa7b_remission_ic = list(
    PNU = c("B-L0"),
    merged = c("B-L0S0")
  ),
  sa8_delta_narrow = list(
    merged = c("B-L0S0")
  ),
  sa8_delta_wide = list(
    merged = c("B-L0S0")
  ),
  sa9_anchor_alt = list(
    merged = c("B-L0S0")
  ),
  sa10_betaMale_direct = list(
    merged = c("B-L0S0")
  ),
  sa11_int_tight = list(
    merged = c("B-L0S0")
  ),
  sa11_int_wide = list(
    merged = c("B-L0S0")
  ),
  sa12_noncure_benchmark = list(
    merged = c("B-L0S0")
  )
)

bundle_settings <- tibble::tribble(
  ~bundle_id,               ~incidence_spec_id,   ~latency_prior_id, ~age_support_trim, ~use_ipcw, ~force_all_susceptible,
  "primary_family",         "I-main",             "L-main",          FALSE,              FALSE,     FALSE,
  "sa1_robust",             "I-robust",           "L-main",          FALSE,              FALSE,     FALSE,
  "sa2_shape_off",          "I-shape-off",        "L-main",          FALSE,              FALSE,     FALSE,
  "sa3_no_external",        "I-noexternal",       "L-main",          FALSE,              FALSE,     FALSE,
  "sa4_age_support",        "I-main",             "L-main",          TRUE,               FALSE,     FALSE,
  "sa5_ltight",             "I-main",             "L-tight",         FALSE,              FALSE,     FALSE,
  "sa5_lwide",              "I-main",             "L-wide",          FALSE,              FALSE,     FALSE,
  "sa7b_remission_ic",      "I-main",             "L-main",          FALSE,              TRUE,      FALSE,
  "sa8_delta_narrow",       "I-delta-narrow",     "L-main",          FALSE,              FALSE,     FALSE,
  "sa8_delta_wide",         "I-delta-wide",       "L-main",          FALSE,              FALSE,     FALSE,
  "sa9_anchor_alt",         "I-anchor-alt",       "L-main",          FALSE,              FALSE,     FALSE,
  "sa10_betaMale_direct",   "I-betaMale-direct",  "L-main",          FALSE,              FALSE,     FALSE,
  "sa11_int_tight",         "I-main",             "L-int-tight",     FALSE,              FALSE,     FALSE,
  "sa11_int_wide",          "I-main",             "L-int-wide",      FALSE,              FALSE,     FALSE,
  "sa12_noncure_benchmark", "I-noexternal",       "L-main",          FALSE,              FALSE,     TRUE
)

# 🔴 Define: helpers for data and priors ===============================
## 🟠 Helpers: I/O wrappers ===============================
safe_write_csv <- function(df, file_path, overwrite = TRUE) {
  if (file.exists(file_path) && !overwrite) stop("파일이 이미 존재합니다: ", file_path)
  readr::write_csv(df, file_path, na = "")
  invisible(file_path)
}

safe_save_rds <- function(object, file_path, overwrite = TRUE) {
  if (file.exists(file_path) && !overwrite) stop("파일이 이미 존재합니다: ", file_path)
  saveRDS(object, file_path, compress = "xz")
  invisible(file_path)
}

## 🟠 Helpers: raw data normalization ===============================
normalize_site_value <- function(x) {
  out <- toupper(trimws(as.character(x)))
  out <- ifelse(grepl("PNU", out), "PNU", out)
  out <- ifelse(grepl("SNU", out), "SNU", out)
  out
}

coalesce_sex_value <- function(sex_num, sex_fact) {
  sx <- suppressWarnings(as.integer(sex_num))
  sf <- toupper(trimws(as.character(sex_fact)))
  sx[is.na(sx) & sf == "FEMALE"] <- 0L
  sx[is.na(sx) & sf == "MALE"] <- 1L
  sx
}

coalesce_age_value <- function(age_exact_entry, age_int) {
  age_exact_num <- suppressWarnings(as.numeric(age_exact_entry))
  age_int_num <- suppressWarnings(as.numeric(age_int))
  out <- age_exact_num
  out[is.na(out)] <- age_int_num[is.na(out)]
  out
}

make_age_group <- function(age_base) {
  dplyr::case_when(
    is.na(age_base) ~ NA_character_,
    age_base < 20 ~ "<20",
    age_base < 30 ~ "20-29",
    age_base >= 30 ~ "30+",
    TRUE ~ NA_character_
  )
}

assert_required_columns <- function(df, cols) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0L) {
    stop("다음 필수 컬럼이 없습니다: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

prepare_master_data <- function(raw_df) {
  required_cols <- c(
    "id", "site", "sex_num", "sex_fact", "age_int", "age_exact_entry",
    "days_followup", "status_num"
  )
  assert_required_columns(raw_df, required_cols)
  
  out <- raw_df %>%
    dplyr::mutate(
      site = normalize_site_value(site),
      sex = coalesce_sex_value(sex_num, sex_fact),
      age_base = coalesce_age_value(age_exact_entry, age_int),
      days_followup = suppressWarnings(as.numeric(days_followup)),
      status_num = suppressWarnings(as.integer(status_num)),
      subject_id = paste0(site, "::", as.character(id)),
      event_primary = dplyr::if_else(status_num == 1L, 1L, 0L, missing = NA_integer_),
      time_days_raw = days_followup,
      time_years = pmax(days_followup / 365.25, 1e-8),
      age_group = make_age_group(age_base),
      sex_label = dplyr::case_when(
        sex == 0L ~ "Female",
        sex == 1L ~ "Male",
        TRUE ~ NA_character_
      )
    )
  
  if (anyDuplicated(out$subject_id) > 0L) {
    dup_ids <- unique(out$subject_id[duplicated(out$subject_id)])
    stop("site + id 중복이 존재합니다: ", paste(dup_ids, collapse = ", "))
  }
  if (any(is.na(out$site))) stop("site 결측이 존재합니다.")
  if (any(!out$site %in% c("PNU", "SNU"))) {
    bad <- unique(out$site[!out$site %in% c("PNU", "SNU")])
    stop("허용되지 않는 site 값이 존재합니다: ", paste(bad, collapse = ", "))
  }
  if (any(is.na(out$status_num))) stop("status_num 결측이 존재합니다.")
  if (any(!out$status_num %in% c(0L, 1L, 2L))) stop("status_num은 0/1/2만 허용됩니다.")
  if (any(is.na(out$days_followup))) stop("days_followup 결측이 존재합니다.")
  if (any(out$days_followup < 0)) stop("days_followup < 0 값이 존재합니다.")
  if (any(out$site == "SNU" & out$status_num == 2L)) stop("SNU subset에 remission(status_num == 2)이 존재합니다.")
  out
}

prepare_branch_data <- function(master_df, cohort_id, age_support_trim = FALSE) {
  df <- switch(
    cohort_id,
    "merged" = master_df,
    "PNU" = dplyr::filter(master_df, site == "PNU"),
    "SNU" = dplyr::filter(master_df, site == "SNU"),
    stop("지원하지 않는 cohort_id: ", cohort_id)
  )
  
  if (age_support_trim) {
    df <- dplyr::filter(df, age_base < 40)
  }
  
  df <- df %>%
    dplyr::filter(
      !is.na(subject_id),
      !is.na(site),
      !is.na(sex),
      !is.na(age_base),
      !is.na(age_group),
      !is.na(time_days_raw),
      !is.na(time_years),
      !is.na(event_primary)
    )
  
  if (nrow(df) == 0L) {
    stop("분석 가능한 행이 없습니다: cohort=", cohort_id, ", age_support_trim=", age_support_trim)
  }
  
  age_center <- mean(df$age_base)
  age_scale_2sd <- 2 * stats::sd(df$age_base)
  if (!is.finite(age_scale_2sd) || age_scale_2sd <= 0) age_scale_2sd <- 1
  
  df <- df %>%
    dplyr::mutate(
      age_s = (age_base - age_center) / age_scale_2sd,
      site_num = dplyr::if_else(site == "SNU", 1L, 0L)
    )
  
  meta <- tibble::tibble(
    cohort_id = cohort_id,
    n = nrow(df),
    age_support_trim = age_support_trim,
    age_center = age_center,
    age_scale_2sd = age_scale_2sd,
    n_transition = sum(df$event_primary == 1L),
    n_censored_primary = sum(df$event_primary == 0L),
    n_remission_original = sum(df$status_num == 2L),
    last_event_year = ifelse(any(df$event_primary == 1L), max(df$time_years[df$event_primary == 1L]), NA_real_),
    max_followup_year = max(df$time_years)
  )
  
  list(data = df, meta = meta)
}

build_latency_matrix <- function(df, cohort_id, latency_spec_id) {
  if (cohort_id %in% c("PNU", "SNU")) {
    if (latency_spec_id == "B-L0") {
      X <- cbind(age_s = df$age_s, sex = df$sex)
    } else if (latency_spec_id == "B-L1") {
      X <- cbind(age_s = df$age_s, sex = df$sex, age_s_sex = df$age_s * df$sex)
    } else {
      stop("PNU/SNU에서 허용되지 않는 latency_spec_id: ", latency_spec_id)
    }
  } else if (cohort_id == "merged") {
    if (latency_spec_id == "B-L0S0") {
      X <- cbind(age_s = df$age_s, sex = df$sex)
    } else if (latency_spec_id == "B-L1S0") {
      X <- cbind(age_s = df$age_s, sex = df$sex, age_s_sex = df$age_s * df$sex)
    } else if (latency_spec_id == "B-L0S1") {
      X <- cbind(age_s = df$age_s, sex = df$sex, site = df$site_num)
    } else if (latency_spec_id == "B-L1S1") {
      X <- cbind(age_s = df$age_s, sex = df$sex, age_s_sex = df$age_s * df$sex, site = df$site_num)
    } else {
      stop("merged에서 허용되지 않는 latency_spec_id: ", latency_spec_id)
    }
  } else {
    stop("지원하지 않는 cohort_id: ", cohort_id)
  }
  
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  list(X = X, coef_names = colnames(X))
}

## 🟠 Helpers: external anchors and bundle settings ===============================
compute_external_info <- function() {
  to_person_year <- function(x) x / 10000
  to_prob <- function(rate_py) 1 - exp(-rate_py)
  to_logit <- function(p) qlogis(p)
  s0_from_ci <- function(lcl, ucl) (log(ucl) - log(lcl)) / (2 * 1.96)
  
  anchor_tbl <- tibble::tribble(
    ~sex,      ~age_band, ~rate_10k, ~adj_irr_lcl, ~adj_irr_ucl,
    "Female",  "<20",      0.69,      NA_real_,     NA_real_,
    "Female",  "20-29",    1.71,      1.66,         3.28,
    "Female",  "30-39",    1.24,      1.55,         3.28,
    "Female",  "40-49",    0.94,      1.19,         2.77,
    "Male",    "<20",      1.05,      NA_real_,     NA_real_,
    "Male",    "20-29",    4.15,      1.74,         3.92,
    "Male",    "30-39",    1.96,      0.94,         2.36,
    "Male",    "40-49",    0.98,      0.47,         1.41
  ) %>%
    dplyr::mutate(
      rate_py = to_person_year(rate_10k),
      p1 = to_prob(rate_py),
      y = to_logit(p1)
    )
  
  y_f_lt20 <- anchor_tbl %>% dplyr::filter(sex == "Female", age_band == "<20") %>% dplyr::pull(y)
  y_f_20   <- anchor_tbl %>% dplyr::filter(sex == "Female", age_band == "20-29") %>% dplyr::pull(y)
  y_f_30   <- anchor_tbl %>% dplyr::filter(sex == "Female", age_band == "30-39") %>% dplyr::pull(y)
  y_f_40   <- anchor_tbl %>% dplyr::filter(sex == "Female", age_band == "40-49") %>% dplyr::pull(y)
  y_m_lt20 <- anchor_tbl %>% dplyr::filter(sex == "Male", age_band == "<20") %>% dplyr::pull(y)
  y_m_20   <- anchor_tbl %>% dplyr::filter(sex == "Male", age_band == "20-29") %>% dplyr::pull(y)
  y_m_30   <- anchor_tbl %>% dplyr::filter(sex == "Male", age_band == "30-39") %>% dplyr::pull(y)
  y_m_40   <- anchor_tbl %>% dplyr::filter(sex == "Male", age_band == "40-49") %>% dplyr::pull(y)
  
  mu_main <- c(
    y_m_lt20 - y_f_lt20,
    y_f_20 - y_f_lt20,
    y_f_30 - y_f_lt20,
    (y_m_20 - y_f_20) - (y_m_lt20 - y_f_lt20),
    (y_m_30 - y_f_30) - (y_m_lt20 - y_f_lt20)
  )
  names(mu_main) <- c("betaMale", "beta20_29", "beta30p", "betaMale20_29", "betaMale30p")
  
  y_f_30_alt <- mean(c(y_f_30, y_f_40))
  y_m_30_alt <- mean(c(y_m_30, y_m_40))
  mu_anchor_alt <- c(
    y_m_lt20 - y_f_lt20,
    y_f_20 - y_f_lt20,
    y_f_30_alt - y_f_lt20,
    (y_m_20 - y_f_20) - (y_m_lt20 - y_f_lt20),
    (y_m_30_alt - y_f_30_alt) - (y_m_lt20 - y_f_lt20)
  )
  names(mu_anchor_alt) <- names(mu_main)
  
  s0_main <- c(
    betaMale = max(
      s0_from_ci(1.66, 3.28),
      s0_from_ci(1.55, 3.28),
      s0_from_ci(1.74, 3.92),
      s0_from_ci(0.94, 2.36)
    ),
    beta20_29 = s0_from_ci(1.66, 3.28),
    beta30p = s0_from_ci(1.55, 3.28),
    betaMale20_29 = sqrt(s0_from_ci(1.74, 3.92)^2 + s0_from_ci(1.66, 3.28)^2),
    betaMale30p = sqrt(s0_from_ci(0.94, 2.36)^2 + s0_from_ci(1.55, 3.28)^2)
  )
  
  s0_betaMale_direct <- s0_from_ci(0.41, 0.69)
  
  list(
    alpha0_gp = y_f_lt20,
    anchor_table = anchor_tbl,
    mu_main = mu_main,
    mu_anchor_alt = mu_anchor_alt,
    s0_main = s0_main,
    s0_betaMale_direct = s0_betaMale_direct
  )
}

external_info <- compute_external_info()

build_prior_settings <- function(bundle_id, external_info) {
  mu_beta <- external_info$mu_main
  s0_beta <- external_info$s0_main
  alpha0_gp <- external_info$alpha0_gp
  incidence_regime_int <- 1L
  delta_sd <- 3.5
  robust_weight <- 0.8
  force_all_susceptible <- FALSE
  
  if (bundle_id == "sa1_robust") {
    incidence_regime_int <- 2L
  } else if (bundle_id == "sa2_shape_off") {
    incidence_regime_int <- 3L
  } else if (bundle_id == "sa3_no_external") {
    incidence_regime_int <- 4L
  } else if (bundle_id == "sa8_delta_narrow") {
    incidence_regime_int <- 1L
    delta_sd <- 2.5
  } else if (bundle_id == "sa8_delta_wide") {
    incidence_regime_int <- 1L
    delta_sd <- 5.0
  } else if (bundle_id == "sa9_anchor_alt") {
    incidence_regime_int <- 1L
    mu_beta <- external_info$mu_anchor_alt
  } else if (bundle_id == "sa10_betaMale_direct") {
    incidence_regime_int <- 1L
    s0_beta["betaMale"] <- external_info$s0_betaMale_direct
  } else if (bundle_id == "sa12_noncure_benchmark") {
    incidence_regime_int <- 4L
    force_all_susceptible <- TRUE
  }
  
  list(
    alpha0_gp = alpha0_gp,
    mu_beta = unname(mu_beta),
    s0_beta = unname(s0_beta),
    incidence_regime_int = incidence_regime_int,
    delta_sd = delta_sd,
    robust_weight = robust_weight,
    force_all_susceptible = force_all_susceptible
  )
}

build_latency_prior_settings <- function(latency_prior_id) {
  switch(
    latency_prior_id,
    "L-main" = list(gamma_sd = 1.0, gamma0_sd = 2.0, log_sigma_mean = log(0.8), log_sigma_sd = 0.35),
    "L-tight" = list(gamma_sd = 0.5, gamma0_sd = 2.0, log_sigma_mean = log(0.8), log_sigma_sd = 0.25),
    "L-wide" = list(gamma_sd = 1.5, gamma0_sd = 2.0, log_sigma_mean = log(0.8), log_sigma_sd = 0.50),
    "L-int-tight" = list(gamma_sd = 1.0, gamma0_sd = 1.5, log_sigma_mean = log(0.8), log_sigma_sd = 0.35),
    "L-int-wide" = list(gamma_sd = 1.0, gamma0_sd = 2.5, log_sigma_mean = log(0.8), log_sigma_sd = 0.35),
    stop("지원하지 않는 latency_prior_id: ", latency_prior_id)
  )
}

build_model_registry <- function(active_bundle_ids, bundle_latency_map, bundle_settings, family_lookup) {
  regs <- list()
  
  for (bundle_id in active_bundle_ids) {
    if (!bundle_id %in% names(bundle_latency_map)) next
    bundle_row <- bundle_settings %>% dplyr::filter(bundle_id == !!bundle_id)
    if (nrow(bundle_row) != 1L) stop("bundle_settings 정의 오류: ", bundle_id)
    latency_map <- bundle_latency_map[[bundle_id]]
    
    for (cohort_id in names(latency_map)) {
      for (latency_spec_id in latency_map[[cohort_id]]) {
        for (i in seq_len(nrow(family_lookup))) {
          regs[[length(regs) + 1L]] <- tibble::tibble(
            bundle_id = bundle_id,
            cohort_id = cohort_id,
            latency_spec_id = latency_spec_id,
            family_id = family_lookup$family_id[i],
            family_name = family_lookup$family_name[i],
            incidence_spec_id = bundle_row$incidence_spec_id,
            latency_prior_id = bundle_row$latency_prior_id,
            age_support_trim = bundle_row$age_support_trim,
            use_ipcw = bundle_row$use_ipcw,
            force_all_susceptible = bundle_row$force_all_susceptible
          )
        }
      }
    }
  }
  
  dplyr::bind_rows(regs) %>%
    dplyr::mutate(
      model_id = paste(bundle_id, cohort_id, latency_spec_id, family_name, sep = "__"),
      model_class = dplyr::if_else(force_all_susceptible, "bayes_nocure_benchmark", "bayes_cure"),
      lane = "bayes_AFT"
    ) %>%
    dplyr::select(
      model_id, bundle_id, cohort_id, latency_spec_id, family_id, family_name,
      incidence_spec_id, latency_prior_id, age_support_trim, use_ipcw,
      force_all_susceptible, model_class, lane
    )
}

model_registry <- build_model_registry(active_bundle_ids, bundle_latency_map, bundle_settings, family_lookup)

## 🟠 Helpers: survival utilities ===============================
step_surv_from_survfit <- function(sf, times) {
  s <- summary(sf, times = times, extend = TRUE)
  tibble::tibble(
    time = s$time,
    surv = s$surv,
    lower = s$lower,
    upper = s$upper,
    n_risk = s$n.risk
  )
}

compute_km_yearly_summary <- function(master_df, years) {
  cohorts <- c("merged", "PNU", "SNU")
  purrr::map_dfr(cohorts, function(cohort_id) {
    prep <- prepare_branch_data(master_df, cohort_id, age_support_trim = FALSE)
    df <- prep$data
    sf <- survival::survfit(survival::Surv(time_years, event_primary) ~ 1, data = df)
    km <- step_surv_from_survfit(sf, years)
    tibble::tibble(
      cohort_id = cohort_id,
      year = years,
      S_KM = km$surv,
      Risk_KM = 1 - km$surv,
      LCL = km$lower,
      UCL = km$upper,
      n_risk = km$n_risk,
      last_event_year = prep$meta$last_event_year,
      max_followup_year = prep$meta$max_followup_year,
      tail_flag = dplyr::if_else(years > prep$meta$last_event_year, TRUE, FALSE, missing = TRUE),
      is_extrapolated = dplyr::if_else(years > prep$meta$max_followup_year, TRUE, FALSE, missing = TRUE)
    )
  })
}

safe_timeROC_auc <- function(time, event, score, t_eval) {
  if (length(unique(stats::na.omit(score))) < 2L) return(NA_real_)
  if (sum(event == 1 & time <= t_eval, na.rm = TRUE) < 2L) return(NA_real_)
  if (sum(time > t_eval, na.rm = TRUE) < 2L) return(NA_real_)
  out <- tryCatch(
    timeROC::timeROC(
      T = time,
      delta = event,
      marker = score,
      cause = 1,
      weighting = "marginal",
      times = t_eval,
      iid = FALSE,
      ROC = FALSE
    ),
    error = function(e) NULL
  )
  if (is.null(out)) return(NA_real_)
  as.numeric(out$AUC[1])
}

compute_ipcw_brier <- function(time, event, surv_prob, t_eval) {
  sf_g <- survival::survfit(survival::Surv(time, 1 - event) ~ 1)
  g_t <- summary(sf_g, times = t_eval, extend = TRUE)$surv
  g_t <- pmax(g_t, 1e-8)
  eps <- 1e-8
  g_tminus <- summary(sf_g, times = pmax(time - eps, 0), extend = TRUE)$surv
  g_tminus <- pmax(g_tminus, 1e-8)
  
  part1 <- as.numeric(time <= t_eval & event == 1L) * ((0 - surv_prob)^2 / g_tminus)
  part2 <- as.numeric(time > t_eval) * ((1 - surv_prob)^2 / g_t)
  mean(part1 + part2, na.rm = TRUE)
}

compute_ibs <- function(year_grid, brier_values) {
  keep <- which(is.finite(year_grid) & is.finite(brier_values))
  if (length(keep) < 2L) return(NA_real_)
  x <- year_grid[keep]
  y <- brier_values[keep]
  dx <- diff(x)
  avg_y <- (head(y, -1L) + tail(y, -1L)) / 2
  sum(dx * avg_y) / (max(x) - min(x))
}

build_calibration_data <- function(time, event, risk, eval_time, n_bins = 10L) {
  if (all(!is.finite(risk))) return(NULL)
  if (length(unique(stats::na.omit(risk))) < 2L) return(NULL)
  bins <- min(n_bins, length(unique(stats::na.omit(risk))))
  if (bins < 2L) return(NULL)
  
  df <- tibble::tibble(time = time, event = event, risk = risk) %>%
    dplyr::mutate(bin = dplyr::ntile(risk, bins))
  
  df %>%
    dplyr::group_by(bin) %>%
    dplyr::group_modify(function(dat, key) {
      sf <- survival::survfit(survival::Surv(time, event) ~ 1, data = dat)
      s_eval <- summary(sf, times = eval_time, extend = TRUE)$surv[1]
      tibble::tibble(
        bin = key$bin,
        n = nrow(dat),
        mean_pred_risk = mean(dat$risk, na.rm = TRUE),
        observed_risk_km = 1 - s_eval
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(year = eval_time)
}

compute_ppc_summary <- function(mean_yearly_df, km_df) {
  joined <- mean_yearly_df %>%
    dplyr::left_join(
      km_df %>% dplyr::select(year, Risk_KM),
      by = "year"
    ) %>%
    dplyr::mutate(abs_diff = abs(meanRisk_mean - Risk_KM))
  
  tibble::tibble(
    ppc_rmse = sqrt(mean((joined$meanRisk_mean - joined$Risk_KM)^2, na.rm = TRUE)),
    tail_fit_absdiff_10 = joined$abs_diff[joined$year == max(joined$year)][1] %||% NA_real_,
    tail_fit_flag = dplyr::if_else((joined$abs_diff[joined$year == max(joined$year)][1] %||% Inf) <= tail_fit_absdiff_threshold, TRUE, FALSE, missing = FALSE)
  )
}

## 🟠 Helpers: IPCW for remission sensitivity ===============================
step_cumhaz_at_times <- function(basehaz_df, times) {
  idx <- findInterval(times, basehaz_df$time)
  out <- numeric(length(times))
  keep <- idx > 0L
  out[keep] <- basehaz_df$hazard[idx[keep]]
  out
}

compute_remission_ipcw_weights <- function(df) {
  if (sum(df$status_num == 2L, na.rm = TRUE) < 2L) {
    return(list(
      weights = rep(1, nrow(df)),
      weight_summary = tibble::tibble(
        mean_weight = 1, sd_weight = 0, min_weight = 1, max_weight = 1,
        p01_weight = 1, p99_weight = 1, remission_events = sum(df$status_num == 2L)
      ),
      censor_model = NULL
    ))
  }
  
  remission_event <- as.integer(df$status_num == 2L)
  
  cox_fit <- survival::coxph(
    survival::Surv(time_years, remission_event) ~ age_s + sex,
    data = df,
    ties = "breslow",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  
  bh <- survival::basehaz(cox_fit, centered = FALSE)
  lp <- stats::predict(cox_fit, type = "lp")
  H0_t <- step_cumhaz_at_times(bh, df$time_years)
  G_i <- exp(-H0_t * exp(lp))
  G_i <- pmax(G_i, 1e-8)
  
  sf_marg <- survival::survfit(survival::Surv(time_years, remission_event) ~ 1, data = df)
  G_marg <- summary(sf_marg, times = df$time_years, extend = TRUE)$surv
  G_marg <- pmax(G_marg, 1e-8)
  
  sw <- G_marg / G_i
  q <- stats::quantile(sw, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  sw <- pmin(pmax(sw, q[1]), q[2])
  
  list(
    weights = sw,
    weight_summary = tibble::tibble(
      mean_weight = mean(sw, na.rm = TRUE),
      sd_weight = stats::sd(sw, na.rm = TRUE),
      min_weight = min(sw, na.rm = TRUE),
      max_weight = max(sw, na.rm = TRUE),
      p01_weight = q[1],
      p99_weight = q[2],
      remission_events = sum(remission_event, na.rm = TRUE)
    ),
    censor_model = cox_fit
  )
}

## 🟠 Helpers: Stan data and diagnostics ===============================
build_stan_data <- function(df, latency_obj, prior_settings, latency_prior_settings, family_id, obs_weight) {
  list(
    N = nrow(df),
    event = as.integer(df$event_primary),
    time = as.vector(df$time_years),
    z = as.integer(df$sex),
    x20 = as.integer(df$age_group == "20-29"),
    x30 = as.integer(df$age_group == "30+"),
    force_all_susceptible = as.integer(prior_settings$force_all_susceptible),
    P = ncol(latency_obj$X),
    X_lat = latency_obj$X,
    family_id = as.integer(family_id),
    incidence_regime = as.integer(prior_settings$incidence_regime_int),
    alpha0_gp = as.numeric(prior_settings$alpha0_gp),
    mu_beta = as.vector(prior_settings$mu_beta),
    s0_beta = as.vector(prior_settings$s0_beta),
    delta_sd = as.numeric(prior_settings$delta_sd),
    kappa_sd = 0.37,
    robust_weight = as.numeric(prior_settings$robust_weight),
    gamma_sd = as.numeric(latency_prior_settings$gamma_sd),
    gamma0_sd = as.numeric(latency_prior_settings$gamma0_sd),
    log_sigma_mean = as.numeric(latency_prior_settings$log_sigma_mean),
    log_sigma_sd = as.numeric(latency_prior_settings$log_sigma_sd),
    obs_weight = as.vector(obs_weight)
  )
}

compute_chain_bfmi <- function(fit) {
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  bfmi_vals <- vapply(sp, function(mat) {
    e <- mat[, "energy__"]
    if (length(e) < 2L || stats::var(e) <= 0) return(NA_real_)
    mean(diff(e)^2) / stats::var(e)
  }, numeric(1))
  bfmi_vals
}

summarize_fit_diagnostics <- function(fit, max_treedepth_used, ess_threshold = 400) {
  draw_summ <- tryCatch(
    posterior::summarise_draws(
      posterior::as_draws_array(fit),
      rhat = posterior::rhat,
      ess_bulk = posterior::ess_bulk,
      ess_tail = posterior::ess_tail
    ),
    error = function(e) NULL
  )
  
  if (is.null(draw_summ) || !all(c("rhat", "ess_bulk", "ess_tail") %in% names(draw_summ))) {
    stan_sum <- rstan::summary(fit)$summary
    max_rhat <- if ("Rhat" %in% colnames(stan_sum)) max(stan_sum[, "Rhat"], na.rm = TRUE) else NA_real_
    min_ess_bulk <- if ("n_eff" %in% colnames(stan_sum)) min(stan_sum[, "n_eff"], na.rm = TRUE) else NA_real_
    min_ess_tail <- min_ess_bulk
  } else {
    max_rhat <- max(draw_summ$rhat, na.rm = TRUE)
    min_ess_bulk <- min(draw_summ$ess_bulk, na.rm = TRUE)
    min_ess_tail <- min(draw_summ$ess_tail, na.rm = TRUE)
  }
  
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergence_count <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
  treedepth_hits <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= max_treedepth_used), numeric(1)))
  bfmi_vals <- compute_chain_bfmi(fit)
  min_bfmi <- suppressWarnings(min(bfmi_vals, na.rm = TRUE))
  if (!is.finite(min_bfmi)) min_bfmi <- NA_real_
  
  needs_refit <- (
    (!is.na(max_rhat) && max_rhat >= 1.01) ||
      (!is.na(min_ess_bulk) && min_ess_bulk < ess_threshold) ||
      (!is.na(min_ess_tail) && min_ess_tail < ess_threshold) ||
      divergence_count > 0L ||
      treedepth_hits > 0L ||
      (!is.na(min_bfmi) && min_bfmi < 0.3)
  )
  
  tibble::tibble(
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    divergence_count = divergence_count,
    treedepth_hits = treedepth_hits,
    min_bfmi = min_bfmi,
    needs_refit = needs_refit
  )
}

sample_with_retry <- function(stan_model_obj, stan_data, seed) {
  fit_try_1 <- tryCatch(
    rstan::sampling(
      object = stan_model_obj,
      data = stan_data,
      seed = seed,
      chains = chains_main,
      iter = warmup_main + iter_sampling_main,
      warmup = warmup_main,
      control = list(adapt_delta = adapt_delta_main, max_treedepth = max_treedepth_main),
      refresh = 100
    ),
    error = function(e) e
  )
  
  if (inherits(fit_try_1, "error")) {
    fit_try_2 <- tryCatch(
      rstan::sampling(
        object = stan_model_obj,
        data = stan_data,
        seed = seed + 1L,
        chains = chains_refit,
        iter = warmup_refit + iter_sampling_refit,
        warmup = warmup_refit,
        control = list(adapt_delta = adapt_delta_refit, max_treedepth = max_treedepth_refit),
        refresh = 100
      ),
      error = function(e) e
    )
    if (inherits(fit_try_2, "error")) {
      return(list(fit = NULL, diagnostics = NULL, refit_used = FALSE, fit_error = paste(fit_try_1$message, fit_try_2$message, sep = " | ")))
    }
    diag2 <- summarize_fit_diagnostics(fit_try_2, max_treedepth_refit, ess_threshold)
    return(list(fit = fit_try_2, diagnostics = diag2, refit_used = TRUE, fit_error = NA_character_))
  }
  
  diag1 <- summarize_fit_diagnostics(fit_try_1, max_treedepth_main, ess_threshold)
  if (!isTRUE(diag1$needs_refit[1])) {
    return(list(fit = fit_try_1, diagnostics = diag1, refit_used = FALSE, fit_error = NA_character_))
  }
  
  fit_try_2 <- tryCatch(
    rstan::sampling(
      object = stan_model_obj,
      data = stan_data,
      seed = seed + 1L,
      chains = chains_refit,
      iter = warmup_refit + iter_sampling_refit,
      warmup = warmup_refit,
      control = list(adapt_delta = adapt_delta_refit, max_treedepth = max_treedepth_refit),
      refresh = 100
    ),
    error = function(e) e
  )
  
  if (inherits(fit_try_2, "error")) {
    return(list(fit = fit_try_1, diagnostics = diag1, refit_used = FALSE, fit_error = paste("refit_failed:", fit_try_2$message)))
  }
  
  diag2 <- summarize_fit_diagnostics(fit_try_2, max_treedepth_refit, ess_threshold)
  list(fit = fit_try_2, diagnostics = diag2, refit_used = TRUE, fit_error = NA_character_)
}

## 🟠 Helpers: posterior prediction ===============================
thin_draw_matrix <- function(draw_mat, max_draws = 1000L) {
  if (nrow(draw_mat) <= max_draws) return(draw_mat)
  idx <- unique(round(seq(1, nrow(draw_mat), length.out = max_draws)))
  draw_mat[idx, , drop = FALSE]
}

extract_core_draws <- function(fit, P, max_draws = prediction_draws_max) {
  pars <- c(
    "alpha0",
    paste0("beta[", 1:5, "]"),
    "log_kappa",
    "gamma0",
    paste0("gamma[", seq_len(P), "]"),
    "log_sigma",
    "delta_equiv"
  )
  draw_mat <- as.matrix(fit, pars = pars)
  thin_draw_matrix(draw_mat, max_draws = max_draws)
}

surv_u_matrix <- function(mu_lat, sigma_vec, year, family_id) {
  S <- nrow(mu_lat)
  N <- ncol(mu_lat)
  
  if (family_id == 1L) {
    return(exp(-year / exp(mu_lat)))
  }
  
  logt_minus_mu <- matrix(log(year), nrow = S, ncol = N) - mu_lat
  scaled <- sweep(logt_minus_mu, 1, sigma_vec, "/")
  
  if (family_id == 2L) {
    z <- exp(scaled)
    return(exp(-z))
  }
  if (family_id == 3L) {
    z <- exp(scaled)
    return(1 / (1 + z))
  }
  if (family_id == 4L) {
    return(pnorm(scaled, lower.tail = FALSE))
  }
  stop("지원하지 않는 family_id: ", family_id)
}

summarize_posterior_vector <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    median = stats::median(x, na.rm = TRUE),
    lcl = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    ucl = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  )
}

predict_from_fit <- function(
    fit,
    df,
    latency_obj,
    family_id,
    years,
    force_all_susceptible = FALSE,
    alpha0_gp = external_info$alpha0_gp
) {
  draw_mat <- extract_core_draws(fit, P = ncol(latency_obj$X), max_draws = prediction_draws_max)
  S <- nrow(draw_mat)
  N <- nrow(df)
  
  alpha0 <- draw_mat[, "alpha0"]
  beta_mat <- draw_mat[, paste0("beta[", 1:5, "]"), drop = FALSE]
  gamma0 <- draw_mat[, "gamma0"]
  gamma_cols <- paste0("gamma[", seq_len(ncol(latency_obj$X)), "]")
  gamma_mat <- draw_mat[, gamma_cols, drop = FALSE]
  log_sigma <- draw_mat[, "log_sigma"]
  log_kappa <- draw_mat[, "log_kappa"]
  
  eta_inc <- matrix(alpha0, nrow = S, ncol = N) +
    tcrossprod(beta_mat[, 1], df$sex) +
    tcrossprod(beta_mat[, 2], as.integer(df$age_group == "20-29")) +
    tcrossprod(beta_mat[, 3], as.integer(df$age_group == "30+")) +
    tcrossprod(beta_mat[, 4], as.integer(df$sex * (df$age_group == "20-29"))) +
    tcrossprod(beta_mat[, 5], as.integer(df$sex * (df$age_group == "30+")))
  
  pi_mat <- if (force_all_susceptible) {
    matrix(1, nrow = S, ncol = N)
  } else {
    plogis(eta_inc)
  }
  
  mu_lat <- matrix(gamma0, nrow = S, ncol = N) + gamma_mat %*% t(latency_obj$X)
  sigma_vec <- if (family_id == 1L) rep(1, S) else exp(log_sigma)
  
  subject_static <- tibble::tibble(
    subject_id = df$subject_id,
    site = df$site,
    sex = df$sex,
    age_base = df$age_base,
    age_group = df$age_group,
    time_years_observed = df$time_years,
    event_primary = df$event_primary,
    susceptible_prob_mean = colMeans(pi_mat),
    susceptible_prob_median = apply(pi_mat, 2, stats::median),
    susceptible_prob_lcl = apply(pi_mat, 2, stats::quantile, probs = 0.025),
    susceptible_prob_ucl = apply(pi_mat, 2, stats::quantile, probs = 0.975),
    cure_prob_mean = colMeans(1 - pi_mat),
    cure_prob_median = apply(1 - pi_mat, 2, stats::median),
    cure_prob_lcl = apply(1 - pi_mat, 2, stats::quantile, probs = 0.025),
    cure_prob_ucl = apply(1 - pi_mat, 2, stats::quantile, probs = 0.975)
  )
  
  yearly_list <- vector("list", length(years))
  mean_risk_draws <- matrix(NA_real_, nrow = S, ncol = length(years))
  colnames(mean_risk_draws) <- paste0("year_", years)
  
  km_tail_cut <- if (any(df$event_primary == 1L)) max(df$time_years[df$event_primary == 1L]) else NA_real_
  max_followup <- max(df$time_years, na.rm = TRUE)
  
  for (k_idx in seq_along(years)) {
    yr <- years[k_idx]
    Su <- surv_u_matrix(mu_lat, sigma_vec, yr, family_id)
    Spop <- if (force_all_susceptible) Su else (1 - pi_mat) + pi_mat * Su
    Risk <- 1 - Spop
    mean_risk_draws[, k_idx] <- rowMeans(Risk)
    
    yearly_list[[k_idx]] <- tibble::tibble(
      subject_id = df$subject_id,
      year = yr,
      survival_mean = colMeans(Spop),
      survival_median = apply(Spop, 2, stats::median),
      survival_lcl = apply(Spop, 2, stats::quantile, probs = 0.025),
      survival_ucl = apply(Spop, 2, stats::quantile, probs = 0.975),
      risk_mean = colMeans(Risk),
      risk_median = apply(Risk, 2, stats::median),
      risk_lcl = apply(Risk, 2, stats::quantile, probs = 0.025),
      risk_ucl = apply(Risk, 2, stats::quantile, probs = 0.975),
      is_extrapolated = yr > max_followup,
      tail_flag = yr > km_tail_cut
    )
  }
  
  subject_yearly <- dplyr::bind_rows(yearly_list)
  
  mean_yearly <- purrr::map_dfr(seq_along(years), function(k_idx) {
    yr <- years[k_idx]
    stat_surv <- summarize_posterior_vector(1 - mean_risk_draws[, k_idx])
    stat_risk <- summarize_posterior_vector(mean_risk_draws[, k_idx])
    tibble::tibble(
      year = yr,
      meanS_mean = stat_surv["mean"],
      meanS_median = stat_surv["median"],
      meanS_lcl = stat_surv["lcl"],
      meanS_ucl = stat_surv["ucl"],
      meanRisk_mean = stat_risk["mean"],
      meanRisk_median = stat_risk["median"],
      meanRisk_lcl = stat_risk["lcl"],
      meanRisk_ucl = stat_risk["ucl"],
      is_extrapolated = yr > max_followup,
      tail_flag = yr > km_tail_cut
    )
  })
  
  cure_frac_draw <- rowMeans(1 - pi_mat)
  cure_frac_summary <- summarize_posterior_vector(cure_frac_draw)
  delta_summary <- summarize_posterior_vector(alpha0 - alpha0_gp)
  kappa_summary <- summarize_posterior_vector(exp(log_kappa))
  sigma_summary <- summarize_posterior_vector(if (family_id == 1L) rep(1, S) else exp(log_sigma))
  
  list(
    subject_static = subject_static,
    subject_yearly = subject_yearly,
    mean_yearly = mean_yearly,
    cure_frac_summary = cure_frac_summary,
    delta_summary = delta_summary,
    kappa_summary = kappa_summary,
    sigma_summary = sigma_summary,
    mean_risk_draws = mean_risk_draws
  )
}

compute_fit_metrics <- function(df, subject_yearly_df, mean_yearly_df, km_df_cohort, fit) {
  yearly_metrics <- purrr::map_dfr(prediction_years, function(yr) {
    subj_row <- subject_yearly_df %>% dplyr::filter(year == !!yr)
    if (nrow(subj_row) == 0L) {
      return(tibble::tibble(year = yr, AUC = NA_real_, Brier = NA_real_))
    }
    auc_val <- safe_timeROC_auc(
      time = df$time_years,
      event = df$event_primary,
      score = subj_row$risk_mean,
      t_eval = yr
    )
    brier_val <- compute_ipcw_brier(
      time = df$time_years,
      event = df$event_primary,
      surv_prob = subj_row$survival_mean,
      t_eval = yr
    )
    tibble::tibble(year = yr, AUC = auc_val, Brier = brier_val)
  })
  
  ibs <- compute_ibs(yearly_metrics$year, yearly_metrics$Brier)
  
  log_lik_mat <- tryCatch(as.matrix(fit, pars = "log_lik"), error = function(e) NULL)
  loo_obj <- if (is.null(log_lik_mat)) NULL else tryCatch(loo::loo(log_lik_mat), error = function(e) NULL)
  waic_obj <- if (is.null(log_lik_mat)) NULL else tryCatch(loo::waic(log_lik_mat), error = function(e) NULL)
  
  ppc <- compute_ppc_summary(mean_yearly_df, km_df_cohort)
  
  mean_yearly_with_km <- mean_yearly_df %>%
    dplyr::left_join(
      km_df_cohort %>% dplyr::select(year, Risk_KM, S_KM, n_risk),
      by = "year"
    ) %>%
    dplyr::mutate(
      Shortfall_model_minus_KM = meanRisk_mean - Risk_KM,
      Ratio_KM_to_model = Risk_KM / meanRisk_mean,
      MissedPos100 = 100 * Shortfall_model_minus_KM
    )
  
  list(
    yearly_metrics = yearly_metrics,
    ibs = ibs,
    looic = if (!is.null(loo_obj)) loo_obj$estimates["looic", "Estimate"] else NA_real_,
    elpd_loo = if (!is.null(loo_obj)) loo_obj$estimates["elpd_loo", "Estimate"] else NA_real_,
    p_loo = if (!is.null(loo_obj)) loo_obj$estimates["p_loo", "Estimate"] else NA_real_,
    waic = if (!is.null(waic_obj)) waic_obj$estimates["waic", "Estimate"] else NA_real_,
    elpd_waic = if (!is.null(waic_obj)) waic_obj$estimates["elpd_waic", "Estimate"] else NA_real_,
    p_waic = if (!is.null(waic_obj)) waic_obj$estimates["p_waic", "Estimate"] else NA_real_,
    ppc = ppc,
    mean_yearly_with_km = mean_yearly_with_km
  )
}

compute_selected_contrast <- function(draw_mat_1, draw_mat_2, years, direction_label) {
  n_draws <- min(nrow(draw_mat_1), nrow(draw_mat_2))
  if (is.null(n_draws) || !is.finite(n_draws) || n_draws < 10L) {
    return(tibble::tibble(
      year = years,
      contrast_label = direction_label,
      diff_mean = NA_real_,
      diff_median = NA_real_,
      diff_lcl = NA_real_,
      diff_ucl = NA_real_,
      p_diff_gt_0 = NA_real_
    ))
  }
  
  idx1 <- seq_len(n_draws)
  idx2 <- seq_len(n_draws)
  
  purrr::map_dfr(seq_along(years), function(j) {
    diff_draw <- drop(draw_mat_1[idx1, j, drop = FALSE]) - drop(draw_mat_2[idx2, j, drop = FALSE])
    tibble::tibble(
      year = years[j],
      contrast_label = direction_label,
      diff_mean = mean(diff_draw),
      diff_median = stats::median(diff_draw),
      diff_lcl = stats::quantile(diff_draw, 0.025),
      diff_ucl = stats::quantile(diff_draw, 0.975),
      p_diff_gt_0 = mean(diff_draw > 0)
    )
  })
}

## 🟠 Helpers: remission competing-risk summaries ===============================
extract_cif_at_times <- function(cif_obj, failcode, times) {
  nm <- grep(paste0("^", failcode), names(cif_obj), value = TRUE)[1]
  if (is.na(nm) || length(nm) == 0L) {
    return(rep(NA_real_, length(times)))
  }
  obj <- cif_obj[[nm]]
  stats::approx(
    x = c(0, obj$time),
    y = c(0, obj$est),
    xout = times,
    method = "constant",
    yleft = 0,
    rule = 2,
    f = 0
  )$y
}

run_competing_risk_summary <- function(master_df, cohort_id, years) {
  prep <- prepare_branch_data(master_df, cohort_id, age_support_trim = FALSE)
  df <- prep$data %>% dplyr::mutate(fstatus = status_num)
  
  if (!any(df$fstatus == 2L)) {
    return(list(
      yearly = tibble::tibble(
        cohort_id = cohort_id,
        year = years,
        transition_cif = NA_real_,
        remission_cif = NA_real_
      ),
      coef = tibble::tibble(
        cohort_id = cohort_id,
        term = character(0),
        estimate = numeric(0),
        se = numeric(0),
        z = numeric(0),
        p_value = numeric(0)
      )
    ))
  }
  
  cif_obj <- cmprsk::cuminc(ftime = df$time_years, fstatus = df$fstatus, cencode = 0)
  yearly <- tibble::tibble(
    cohort_id = cohort_id,
    year = years,
    transition_cif = extract_cif_at_times(cif_obj, failcode = 1, times = years),
    remission_cif = extract_cif_at_times(cif_obj, failcode = 2, times = years)
  )
  
  X_fg <- as.matrix(df %>% dplyr::transmute(age_s = age_s, sex = sex))
  fg_fit <- tryCatch(
    cmprsk::crr(ftime = df$time_years, fstatus = df$fstatus, cov1 = X_fg, failcode = 1, cencode = 0),
    error = function(e) NULL
  )
  
  coef_tbl <- if (is.null(fg_fit)) {
    tibble::tibble(
      cohort_id = cohort_id,
      term = character(0),
      estimate = numeric(0),
      se = numeric(0),
      z = numeric(0),
      p_value = numeric(0)
    )
  } else {
    se <- sqrt(diag(fg_fit$var))
    z <- fg_fit$coef / se
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    tibble::tibble(
      cohort_id = cohort_id,
      term = names(fg_fit$coef),
      estimate = unname(fg_fit$coef),
      se = unname(se),
      z = unname(z),
      p_value = unname(p)
    )
  }
  
  list(yearly = yearly, coef = coef_tbl)
}

# 🔴 Load: raw input file ===============================
raw_data <- readr::read_csv(data_path, show_col_types = FALSE, progress = FALSE)

# 🔴 Transform: cleaned analysis dataset ===============================
master_data <- prepare_master_data(raw_data)

# 🔴 Summarize: cohort QC and KM benchmark ===============================
branch_qc <- purrr::map_dfr(c("merged", "PNU", "SNU"), function(cohort_id) {
  prepare_branch_data(master_data, cohort_id, age_support_trim = FALSE)$meta
})

km_yearly_summary <- compute_km_yearly_summary(master_data, prediction_years)

# 🔴 Compile: inline Stan model ===============================
stan_code <- "
functions {
  real log_surv_u(real t, real mu, real sigma, int family_id) {
    real z;
    if (family_id == 1) {
      return -t / exp(mu);
    } else if (family_id == 2) {
      z = exp((log(t) - mu) / sigma);
      return -z;
    } else if (family_id == 3) {
      z = exp((log(t) - mu) / sigma);
      return -log1p(z);
    } else if (family_id == 4) {
      return normal_lccdf(log(t) | mu, sigma);
    } else {
      reject(\"family_id must be 1,2,3,4\");
    }
    return negative_infinity();
  }

  real log_pdf_u(real t, real mu, real sigma, int family_id) {
    real z;
    if (family_id == 1) {
      return -mu - t / exp(mu);
    } else if (family_id == 2) {
      z = exp((log(t) - mu) / sigma);
      return log(z) - log(t) - log(sigma) - z;
    } else if (family_id == 3) {
      z = (log(t) - mu) / sigma;
      return logistic_lpdf(z | 0, 1) - log(sigma) - log(t);
    } else if (family_id == 4) {
      return normal_lpdf(log(t) | mu, sigma) - log(t);
    } else {
      reject(\"family_id must be 1,2,3,4\");
    }
    return negative_infinity();
  }

  real cure_loglik_i(real t, int event, real pi, real mu, real sigma, int family_id) {
    if (event == 1) {
      return log(pi) + log_pdf_u(t, mu, sigma, family_id);
    } else {
      return log_sum_exp(log1m(pi), log(pi) + log_surv_u(t, mu, sigma, family_id));
    }
  }
}
data {
  int<lower=1> N;
  array[N] int<lower=0,upper=1> event;
  vector<lower=1e-8>[N] time;
  array[N] int<lower=0,upper=1> z;
  array[N] int<lower=0,upper=1> x20;
  array[N] int<lower=0,upper=1> x30;
  int<lower=0,upper=1> force_all_susceptible;
  int<lower=1> P;
  matrix[N,P] X_lat;
  int<lower=1,upper=4> family_id;
  int<lower=1,upper=4> incidence_regime;
  real alpha0_gp;
  vector[5] mu_beta;
  vector<lower=0>[5] s0_beta;
  real<lower=0> delta_sd;
  real<lower=0> kappa_sd;
  real<lower=0,upper=1> robust_weight;
  real<lower=0> gamma_sd;
  real<lower=0> gamma0_sd;
  real log_sigma_mean;
  real<lower=0> log_sigma_sd;
  vector<lower=0>[N] obs_weight;
}
parameters {
  real alpha0;
  vector[5] beta;
  real log_kappa;
  real gamma0;
  vector[P] gamma;
  real log_sigma;
}
transformed parameters {
  real<lower=0> kappa = exp(log_kappa);
  real<lower=0> sigma;
  vector[N] eta_inc;
  vector[N] mu_lat;
  sigma = family_id == 1 ? 1.0 : exp(log_sigma);

  for (i in 1:N) {
    eta_inc[i] = alpha0
      + beta[1] * z[i]
      + beta[2] * x20[i]
      + beta[3] * x30[i]
      + beta[4] * (z[i] * x20[i])
      + beta[5] * (z[i] * x30[i]);
  }

  mu_lat = gamma0 + X_lat * gamma;
}
model {
  if (incidence_regime == 1 || incidence_regime == 2 || incidence_regime == 3) {
    alpha0 ~ normal(alpha0_gp, delta_sd);
  } else {
    alpha0 ~ normal(0, 2.5);
  }

  if (incidence_regime == 1) {
    log_kappa ~ normal(0, kappa_sd);
    for (j in 1:5) {
      beta[j] ~ normal(mu_beta[j], kappa * s0_beta[j]);
    }
  } else if (incidence_regime == 2) {
    log_kappa ~ normal(0, kappa_sd);
    for (j in 1:5) {
      target += log_mix(
        robust_weight,
        normal_lpdf(beta[j] | mu_beta[j], kappa * s0_beta[j]),
        normal_lpdf(beta[j] | 0, 1)
      );
    }
  } else if (incidence_regime == 3) {
    log_kappa ~ normal(0, 0.001);
    beta ~ normal(0, 1);
  } else if (incidence_regime == 4) {
    log_kappa ~ normal(0, 0.001);
    beta ~ normal(0, 1);
  }

  gamma0 ~ normal(0, gamma0_sd);
  gamma ~ normal(0, gamma_sd);

  if (family_id == 1) {
    log_sigma ~ normal(0, 0.001);
  } else {
    log_sigma ~ normal(log_sigma_mean, log_sigma_sd);
  }

  for (i in 1:N) {
    real ll;
    real pi_i;
    pi_i = inv_logit(eta_inc[i]);

    if (force_all_susceptible == 1) {
      if (event[i] == 1) {
        ll = log_pdf_u(time[i], mu_lat[i], sigma, family_id);
      } else {
        ll = log_surv_u(time[i], mu_lat[i], sigma, family_id);
      }
    } else {
      ll = cure_loglik_i(time[i], event[i], pi_i, mu_lat[i], sigma, family_id);
    }

    target += obs_weight[i] * ll;
  }
}
generated quantities {
  vector[N] log_lik;
  real delta_equiv;
  delta_equiv = alpha0 - alpha0_gp;

  for (i in 1:N) {
    real ll;
    real pi_i;
    pi_i = inv_logit(eta_inc[i]);

    if (force_all_susceptible == 1) {
      if (event[i] == 1) {
        ll = log_pdf_u(time[i], mu_lat[i], sigma, family_id);
      } else {
        ll = log_surv_u(time[i], mu_lat[i], sigma, family_id);
      }
    } else {
      ll = cure_loglik_i(time[i], event[i], pi_i, mu_lat[i], sigma, family_id);
    }

    log_lik[i] = obs_weight[i] * ll;
  }
}
"

compiled_stan_model <- rstan::stan_model(
  model_code = stan_code,
  model_name = paste0("chr_bayes_cure_inline_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)

# 🔴 Fit: Bayesian grid ===============================
fit_store <- list()
registry_results <- list()
subject_static_store <- list()
subject_yearly_store <- list()
mean_yearly_store <- list()
metrics_store <- list()
mcmc_store <- list()
calibration_store <- list()
cohort_draw_store <- list()
weight_summary_store <- list()

subject_static_template <- tibble::tibble(
  model_id = character(),
  bundle_id = character(),
  cohort_id = character(),
  family_name = character(),
  latency_spec_id = character(),
  incidence_spec_id = character(),
  latency_prior_id = character(),
  subject_id = character(),
  site = character(),
  sex = integer(),
  age_base = double(),
  age_group = character(),
  time_years_observed = double(),
  event_primary = integer(),
  susceptible_prob_mean = double(),
  susceptible_prob_median = double(),
  susceptible_prob_lcl = double(),
  susceptible_prob_ucl = double(),
  cure_prob_mean = double(),
  cure_prob_median = double(),
  cure_prob_lcl = double(),
  cure_prob_ucl = double()
)

subject_yearly_template <- tibble::tibble(
  model_id = character(),
  bundle_id = character(),
  cohort_id = character(),
  family_name = character(),
  latency_spec_id = character(),
  incidence_spec_id = character(),
  latency_prior_id = character(),
  subject_id = character(),
  site = character(),
  sex = integer(),
  age_base = double(),
  age_group = character(),
  year = double(),
  survival_mean = double(),
  survival_median = double(),
  survival_lcl = double(),
  survival_ucl = double(),
  risk_mean = double(),
  risk_median = double(),
  risk_lcl = double(),
  risk_ucl = double(),
  is_extrapolated = logical(),
  tail_flag = logical()
)

mean_yearly_template <- tibble::tibble(
  model_id = character(),
  bundle_id = character(),
  cohort_id = character(),
  family_name = character(),
  latency_spec_id = character(),
  incidence_spec_id = character(),
  latency_prior_id = character(),
  year = double(),
  meanS_mean = double(),
  meanS_median = double(),
  meanS_lcl = double(),
  meanS_ucl = double(),
  meanRisk_mean = double(),
  meanRisk_median = double(),
  meanRisk_lcl = double(),
  meanRisk_ucl = double(),
  is_extrapolated = logical(),
  tail_flag = logical(),
  Risk_KM = double(),
  S_KM = double(),
  n_risk = double(),
  Shortfall_model_minus_KM = double(),
  Ratio_KM_to_model = double(),
  MissedPos100 = double()
)

metrics_template <- tibble::tibble(
  model_id = character(),
  bundle_id = character(),
  cohort_id = character(),
  family_name = character(),
  latency_spec_id = character(),
  incidence_spec_id = character(),
  latency_prior_id = character(),
  year = double(),
  AUC = double(),
  Brier = double(),
  IBS = double(),
  looic = double(),
  waic = double()
)

mcmc_template <- tibble::tibble(
  model_id = character(),
  bundle_id = character(),
  cohort_id = character(),
  family_name = character(),
  latency_spec_id = character(),
  incidence_spec_id = character(),
  latency_prior_id = character(),
  max_rhat = double(),
  min_ess_bulk = double(),
  min_ess_tail = double(),
  divergence_count = double(),
  treedepth_hits = double(),
  min_bfmi = double(),
  needs_refit = logical()
)

calibration_template <- tibble::tibble(
  model_id = character(),
  bundle_id = character(),
  cohort_id = character(),
  family_name = character(),
  latency_spec_id = character(),
  incidence_spec_id = character(),
  latency_prior_id = character(),
  bin = integer(),
  n = integer(),
  mean_pred_risk = double(),
  observed_risk_km = double(),
  year = double()
)

for (row_idx in seq_len(nrow(model_registry))) {
  row_info <- model_registry[row_idx, ]
  
  prep <- prepare_branch_data(
    master_df = master_data,
    cohort_id = row_info$cohort_id,
    age_support_trim = row_info$age_support_trim
  )
  
  df <- prep$data
  latency_obj <- build_latency_matrix(df, row_info$cohort_id, row_info$latency_spec_id)
  prior_settings <- build_prior_settings(row_info$bundle_id, external_info)
  latency_prior_settings <- build_latency_prior_settings(row_info$latency_prior_id)
  
  if (isTRUE(row_info$use_ipcw)) {
    ipcw_obj <- compute_remission_ipcw_weights(df)
    obs_weight <- ipcw_obj$weights
    weight_summary_store[[row_info$model_id]] <- ipcw_obj$weight_summary %>%
      dplyr::mutate(model_id = row_info$model_id)
  } else {
    ipcw_obj <- NULL
    obs_weight <- rep(1, nrow(df))
    weight_summary_store[[row_info$model_id]] <- tibble::tibble(
      model_id = row_info$model_id,
      mean_weight = 1,
      sd_weight = 0,
      min_weight = 1,
      max_weight = 1,
      p01_weight = 1,
      p99_weight = 1,
      remission_events = sum(df$status_num == 2L, na.rm = TRUE)
    )
  }
  
  stan_data <- build_stan_data(
    df = df,
    latency_obj = latency_obj,
    prior_settings = prior_settings,
    latency_prior_settings = latency_prior_settings,
    family_id = row_info$family_id,
    obs_weight = obs_weight
  )
  
  fit_result <- sample_with_retry(
    stan_model_obj = compiled_stan_model,
    stan_data = stan_data,
    seed = sampling_seed + row_idx
  )
  
  if (is.null(fit_result$fit)) {
    registry_results[[row_idx]] <- row_info %>%
      dplyr::mutate(
        fit_status = "failed",
        refit_used = FALSE,
        fit_error = fit_result$fit_error,
        n_analysis = nrow(df),
        n_transition = sum(df$event_primary == 1L),
        n_censored_primary = sum(df$event_primary == 0L),
        age_center = prep$meta$age_center,
        age_scale_2sd = prep$meta$age_scale_2sd,
        delta_mean = NA_real_,
        delta_median = NA_real_,
        delta_lcl = NA_real_,
        delta_ucl = NA_real_,
        kappa_mean = NA_real_,
        kappa_median = NA_real_,
        kappa_lcl = NA_real_,
        kappa_ucl = NA_real_,
        sigma_mean = NA_real_,
        sigma_median = NA_real_,
        sigma_lcl = NA_real_,
        sigma_ucl = NA_real_,
        posterior_cure_fraction_mean = NA_real_,
        posterior_cure_fraction_median = NA_real_,
        posterior_cure_fraction_lcl = NA_real_,
        posterior_cure_fraction_ucl = NA_real_,
        looic = NA_real_,
        elpd_loo = NA_real_,
        p_loo = NA_real_,
        waic = NA_real_,
        elpd_waic = NA_real_,
        p_waic = NA_real_,
        IBS = NA_real_,
        ppc_rmse = NA_real_,
        tail_fit_absdiff_10 = NA_real_,
        tail_fit_flag = NA
      )
    next
  }
  
  fit_obj <- fit_result$fit
  fit_diag <- fit_result$diagnostics
  
  fit_store[[row_info$model_id]] <- list(
    fit = if (save_stanfit_inside_master_rds) fit_obj else NULL,
    data_used = df,
    design_info = list(
      cohort_id = row_info$cohort_id,
      latency_spec_id = row_info$latency_spec_id,
      latency_coef_names = latency_obj$coef_names,
      age_center = prep$meta$age_center,
      age_scale_2sd = prep$meta$age_scale_2sd,
      age_support_trim = row_info$age_support_trim
    ),
    prior_settings = prior_settings,
    latency_prior_settings = latency_prior_settings,
    obs_weight = obs_weight,
    ipcw_censor_model = ipcw_obj$censor_model %||% NULL
  )
  
  pred_cache <- predict_from_fit(
    fit = fit_obj,
    df = df,
    latency_obj = latency_obj,
    family_id = row_info$family_id,
    years = prediction_years,
    force_all_susceptible = row_info$force_all_susceptible,
    alpha0_gp = external_info$alpha0_gp
  )
  
  km_df_cohort <- km_yearly_summary %>% dplyr::filter(cohort_id == row_info$cohort_id)
  fit_metrics <- compute_fit_metrics(
    df = df,
    subject_yearly_df = pred_cache$subject_yearly,
    mean_yearly_df = pred_cache$mean_yearly,
    km_df_cohort = km_df_cohort,
    fit = fit_obj
  )
  
  registry_results[[row_idx]] <- row_info %>%
    dplyr::mutate(
      fit_status = "ok",
      refit_used = fit_result$refit_used,
      fit_error = fit_result$fit_error,
      n_analysis = nrow(df),
      n_transition = sum(df$event_primary == 1L),
      n_censored_primary = sum(df$event_primary == 0L),
      age_center = prep$meta$age_center,
      age_scale_2sd = prep$meta$age_scale_2sd,
      delta_mean = pred_cache$delta_summary["mean"],
      delta_median = pred_cache$delta_summary["median"],
      delta_lcl = pred_cache$delta_summary["lcl"],
      delta_ucl = pred_cache$delta_summary["ucl"],
      kappa_mean = pred_cache$kappa_summary["mean"],
      kappa_median = pred_cache$kappa_summary["median"],
      kappa_lcl = pred_cache$kappa_summary["lcl"],
      kappa_ucl = pred_cache$kappa_summary["ucl"],
      sigma_mean = pred_cache$sigma_summary["mean"],
      sigma_median = pred_cache$sigma_summary["median"],
      sigma_lcl = pred_cache$sigma_summary["lcl"],
      sigma_ucl = pred_cache$sigma_summary["ucl"],
      posterior_cure_fraction_mean = pred_cache$cure_frac_summary["mean"],
      posterior_cure_fraction_median = pred_cache$cure_frac_summary["median"],
      posterior_cure_fraction_lcl = pred_cache$cure_frac_summary["lcl"],
      posterior_cure_fraction_ucl = pred_cache$cure_frac_summary["ucl"],
      looic = fit_metrics$looic,
      elpd_loo = fit_metrics$elpd_loo,
      p_loo = fit_metrics$p_loo,
      waic = fit_metrics$waic,
      elpd_waic = fit_metrics$elpd_waic,
      p_waic = fit_metrics$p_waic,
      IBS = fit_metrics$ibs,
      ppc_rmse = fit_metrics$ppc$ppc_rmse,
      tail_fit_absdiff_10 = fit_metrics$ppc$tail_fit_absdiff_10,
      tail_fit_flag = fit_metrics$ppc$tail_fit_flag
    )
  
  subject_static_store[[row_info$model_id]] <- pred_cache$subject_static %>%
    dplyr::mutate(
      model_id = row_info$model_id,
      bundle_id = row_info$bundle_id,
      cohort_id = row_info$cohort_id,
      family_name = row_info$family_name,
      latency_spec_id = row_info$latency_spec_id,
      incidence_spec_id = row_info$incidence_spec_id,
      latency_prior_id = row_info$latency_prior_id
    ) %>%
    dplyr::select(
      model_id, bundle_id, cohort_id, family_name,
      latency_spec_id, incidence_spec_id, latency_prior_id,
      subject_id, site, sex, age_base, age_group, time_years_observed,
      event_primary, susceptible_prob_mean, susceptible_prob_median,
      susceptible_prob_lcl, susceptible_prob_ucl, cure_prob_mean,
      cure_prob_median, cure_prob_lcl, cure_prob_ucl
    )
  
  subject_yearly_store[[row_info$model_id]] <- pred_cache$subject_yearly %>%
    dplyr::left_join(
      pred_cache$subject_static %>% dplyr::select(subject_id, site, sex, age_base, age_group),
      by = "subject_id"
    ) %>%
    dplyr::mutate(
      model_id = row_info$model_id,
      bundle_id = row_info$bundle_id,
      cohort_id = row_info$cohort_id,
      family_name = row_info$family_name,
      latency_spec_id = row_info$latency_spec_id,
      incidence_spec_id = row_info$incidence_spec_id,
      latency_prior_id = row_info$latency_prior_id
    ) %>%
    dplyr::select(
      model_id, bundle_id, cohort_id, family_name,
      latency_spec_id, incidence_spec_id, latency_prior_id,
      subject_id, site, sex, age_base, age_group,
      year, survival_mean, survival_median, survival_lcl, survival_ucl,
      risk_mean, risk_median, risk_lcl, risk_ucl, is_extrapolated, tail_flag
    )
  
  mean_yearly_store[[row_info$model_id]] <- fit_metrics$mean_yearly_with_km %>%
    dplyr::mutate(
      model_id = row_info$model_id,
      bundle_id = row_info$bundle_id,
      cohort_id = row_info$cohort_id,
      family_name = row_info$family_name,
      latency_spec_id = row_info$latency_spec_id,
      incidence_spec_id = row_info$incidence_spec_id,
      latency_prior_id = row_info$latency_prior_id
    ) %>%
    dplyr::select(
      model_id, bundle_id, cohort_id, family_name,
      latency_spec_id, incidence_spec_id, latency_prior_id,
      year, meanS_mean, meanS_median, meanS_lcl, meanS_ucl,
      meanRisk_mean, meanRisk_median, meanRisk_lcl, meanRisk_ucl,
      is_extrapolated, tail_flag, Risk_KM, S_KM, n_risk,
      Shortfall_model_minus_KM, Ratio_KM_to_model, MissedPos100
    )
  
  metrics_store[[row_info$model_id]] <- fit_metrics$yearly_metrics %>%
    dplyr::mutate(
      model_id = row_info$model_id,
      bundle_id = row_info$bundle_id,
      cohort_id = row_info$cohort_id,
      family_name = row_info$family_name,
      latency_spec_id = row_info$latency_spec_id,
      incidence_spec_id = row_info$incidence_spec_id,
      latency_prior_id = row_info$latency_prior_id,
      IBS = fit_metrics$ibs,
      looic = fit_metrics$looic,
      waic = fit_metrics$waic
    ) %>%
    dplyr::select(
      model_id, bundle_id, cohort_id, family_name,
      latency_spec_id, incidence_spec_id, latency_prior_id,
      year, AUC, Brier, IBS, looic, waic
    )
  
  mcmc_store[[row_info$model_id]] <- fit_diag %>%
    dplyr::mutate(
      model_id = row_info$model_id,
      bundle_id = row_info$bundle_id,
      cohort_id = row_info$cohort_id,
      family_name = row_info$family_name,
      latency_spec_id = row_info$latency_spec_id,
      incidence_spec_id = row_info$incidence_spec_id,
      latency_prior_id = row_info$latency_prior_id
    ) %>%
    dplyr::select(
      model_id, bundle_id, cohort_id, family_name,
      latency_spec_id, incidence_spec_id, latency_prior_id,
      max_rhat, min_ess_bulk, min_ess_tail, divergence_count,
      treedepth_hits, min_bfmi, needs_refit
    )
  
  calibration_store[[row_info$model_id]] <- purrr::map_dfr(calibration_years, function(yr) {
    risk_vec <- subject_yearly_store[[row_info$model_id]] %>%
      dplyr::filter(year == !!yr) %>%
      dplyr::pull(risk_mean)
    
    cal <- build_calibration_data(df$time_years, df$event_primary, risk_vec, yr, n_bins = 10L)
    if (is.null(cal)) return(NULL)
    
    cal %>%
      dplyr::mutate(
        model_id = row_info$model_id,
        bundle_id = row_info$bundle_id,
        cohort_id = row_info$cohort_id,
        family_name = row_info$family_name,
        latency_spec_id = row_info$latency_spec_id,
        incidence_spec_id = row_info$incidence_spec_id,
        latency_prior_id = row_info$latency_prior_id
      ) %>%
      dplyr::select(
        model_id, bundle_id, cohort_id, family_name,
        latency_spec_id, incidence_spec_id, latency_prior_id,
        bin, n, mean_pred_risk, observed_risk_km, year
      )
  })
  
  cohort_draw_store[[row_info$model_id]] <- list(
    mean_risk_draws = pred_cache$mean_risk_draws,
    years = prediction_years
  )
  
  gc()
}

model_registry_final <- bind_rows_or_template(registry_results, model_registry %>% dplyr::mutate(
  fit_status = character(),
  refit_used = logical(),
  fit_error = character(),
  n_analysis = integer(),
  n_transition = integer(),
  n_censored_primary = integer(),
  age_center = double(),
  age_scale_2sd = double(),
  delta_mean = double(),
  delta_median = double(),
  delta_lcl = double(),
  delta_ucl = double(),
  kappa_mean = double(),
  kappa_median = double(),
  kappa_lcl = double(),
  kappa_ucl = double(),
  sigma_mean = double(),
  sigma_median = double(),
  sigma_lcl = double(),
  sigma_ucl = double(),
  posterior_cure_fraction_mean = double(),
  posterior_cure_fraction_median = double(),
  posterior_cure_fraction_lcl = double(),
  posterior_cure_fraction_ucl = double(),
  looic = double(),
  elpd_loo = double(),
  p_loo = double(),
  waic = double(),
  elpd_waic = double(),
  p_waic = double(),
  IBS = double(),
  ppc_rmse = double(),
  tail_fit_absdiff_10 = double(),
  tail_fit_flag = logical()
))

bayes_subject_static <- bind_rows_or_template(subject_static_store, subject_static_template)
bayes_subject_predictions <- bind_rows_or_template(subject_yearly_store, subject_yearly_template)
bayes_yearly_predictions <- bind_rows_or_template(mean_yearly_store, mean_yearly_template)
bayes_metrics <- bind_rows_or_template(metrics_store, metrics_template)
bayes_mcmc_diagnostics <- bind_rows_or_template(mcmc_store, mcmc_template)
bayes_calibration <- bind_rows_or_template(calibration_store, calibration_template)
bayes_weight_summary <- bind_rows_or_template(
  weight_summary_store,
  tibble::tibble(
    model_id = character(),
    mean_weight = double(),
    sd_weight = double(),
    min_weight = double(),
    max_weight = double(),
    p01_weight = double(),
    p99_weight = double(),
    remission_events = integer()
  )
)

# 🔴 Run: remission competing-risk summaries ===============================
if (isTRUE(run_remission_competing_risk)) {
  remission_cr_results <- purrr::map(c("merged", "PNU"), ~ run_competing_risk_summary(master_data, .x, prediction_years))
  names(remission_cr_results) <- c("merged", "PNU")
  remission_cr_yearly <- dplyr::bind_rows(purrr::map(remission_cr_results, "yearly"))
  remission_crr_coef <- dplyr::bind_rows(purrr::map(remission_cr_results, "coef"))
} else {
  remission_cr_results <- list()
  remission_cr_yearly <- tibble::tibble()
  remission_crr_coef <- tibble::tibble()
}

# 🔴 Compare: prior utility and posterior contrasts ===============================
## 🟠 Compare: prior utility matrix ===============================
selected_year_summary <- bayes_yearly_predictions %>%
  dplyr::filter(year %in% selected_contrast_years) %>%
  dplyr::mutate(CrI_width_meanRisk = meanRisk_ucl - meanRisk_lcl) %>%
  dplyr::select(
    model_id, year, cohort_id, family_name, latency_spec_id,
    meanS_mean, meanRisk_mean, meanRisk_lcl, meanRisk_ucl,
    Risk_KM, Shortfall_model_minus_KM, Ratio_KM_to_model,
    MissedPos100, CrI_width_meanRisk
  )

utility_base <- model_registry_final %>%
  dplyr::filter(fit_status == "ok") %>%
  dplyr::select(
    model_id, bundle_id, cohort_id, family_name, latency_spec_id,
    looic, waic, ppc_rmse, posterior_cure_fraction_mean,
    delta_mean, kappa_mean
  )

utility_base_selected <- utility_base %>%
  tidyr::crossing(year = selected_contrast_years) %>%
  dplyr::left_join(selected_year_summary, by = c("model_id", "year"))

get_model_key <- function(bundle_id, cohort_id, family_name, latency_spec_id) {
  paste(bundle_id, cohort_id, latency_spec_id, family_name, sep = "__")
}

extract_utility_row <- function(model_id, utility_df, year) {
  utility_df %>% dplyr::filter(model_id == !!model_id, year == !!year)
}

compute_utility_pairs <- function(cohort_id, family_name, latency_spec_id, year, utility_df) {
  main_key <- get_model_key("primary_family", cohort_id, family_name, latency_spec_id)
  robust_key <- get_model_key("sa1_robust", cohort_id, family_name, latency_spec_id)
  shape_key <- get_model_key("sa2_shape_off", cohort_id, family_name, latency_spec_id)
  noext_key <- get_model_key("sa3_no_external", cohort_id, family_name, latency_spec_id)
  
  main_row <- extract_utility_row(main_key, utility_df, year)
  robust_row <- extract_utility_row(robust_key, utility_df, year)
  shape_row <- extract_utility_row(shape_key, utility_df, year)
  noext_row <- extract_utility_row(noext_key, utility_df, year)
  
  tibble::tibble(
    cohort_id = cohort_id,
    family_name = family_name,
    latency_spec_id = latency_spec_id,
    year = year,
    main_model_id = main_key,
    robust_model_id = robust_key,
    shape_model_id = shape_key,
    noexternal_model_id = noext_key,
    delta_looic_main_minus_noext = (main_row$looic[1] %||% NA_real_) - (noext_row$looic[1] %||% NA_real_),
    delta_looic_robust_minus_noext = (robust_row$looic[1] %||% NA_real_) - (noext_row$looic[1] %||% NA_real_),
    delta_risk_main_minus_noext = (main_row$meanRisk_mean[1] %||% NA_real_) - (noext_row$meanRisk_mean[1] %||% NA_real_),
    delta_risk_robust_minus_noext = (robust_row$meanRisk_mean[1] %||% NA_real_) - (noext_row$meanRisk_mean[1] %||% NA_real_),
    delta_risk_shape_minus_main = (shape_row$meanRisk_mean[1] %||% NA_real_) - (main_row$meanRisk_mean[1] %||% NA_real_),
    widthratio_main_to_noext = (main_row$CrI_width_meanRisk[1] %||% NA_real_) / (noext_row$CrI_width_meanRisk[1] %||% NA_real_),
    kappa_main = main_row$kappa_mean[1] %||% NA_real_,
    kappa_robust = robust_row$kappa_mean[1] %||% NA_real_,
    delta_main = main_row$delta_mean[1] %||% NA_real_,
    delta_robust = robust_row$delta_mean[1] %||% NA_real_,
    ppc_rmse_main = main_row$ppc_rmse[1] %||% NA_real_,
    ppc_rmse_noext = noext_row$ppc_rmse[1] %||% NA_real_
  )
}

utility_combos <- tibble::tribble(
  ~cohort_id, ~latency_spec_id,
  "PNU", "B-L0",
  "SNU", "B-L0",
  "merged", "B-L0S0",
  "merged", "B-L0S1"
)

bayes_prior_utility <- purrr::pmap_dfr(
  tidyr::crossing(
    utility_combos,
    family_name = family_lookup$family_name,
    year = selected_contrast_years
  ),
  function(cohort_id, latency_spec_id, family_name, year) {
    compute_utility_pairs(cohort_id, family_name, latency_spec_id, year, utility_base_selected)
  }
)

## 🟠 Compare: posterior contrasts ===============================
contrast_rows <- list()

for (family_name in family_lookup$family_name) {
  for (yr in selected_contrast_years) {
    for (cohort_id in c("PNU", "SNU")) {
      id_l0 <- get_model_key("primary_family", cohort_id, family_name, "B-L0")
      id_l1 <- get_model_key("primary_family", cohort_id, family_name, "B-L1")
      if (!is.null(cohort_draw_store[[id_l0]]) && !is.null(cohort_draw_store[[id_l1]])) {
        col_idx <- match(yr, cohort_draw_store[[id_l1]]$years)
        contrast_rows[[length(contrast_rows) + 1L]] <- compute_selected_contrast(
          draw_mat_1 = cohort_draw_store[[id_l1]]$mean_risk_draws[, col_idx, drop = FALSE],
          draw_mat_2 = cohort_draw_store[[id_l0]]$mean_risk_draws[, col_idx, drop = FALSE],
          years = yr,
          direction_label = paste0("interaction_L1_minus_L0__", cohort_id, "__", family_name)
        ) %>%
          dplyr::mutate(model_id_1 = id_l1, model_id_2 = id_l0)
      }
    }
    
    merged_pairs <- list(
      c("B-L1S0", "B-L0S0", "interaction_L1_minus_L0__merged__S0"),
      c("B-L1S1", "B-L0S1", "interaction_L1_minus_L0__merged__S1"),
      c("B-L0S1", "B-L0S0", "site_S1_minus_S0__merged__L0"),
      c("B-L1S1", "B-L1S0", "site_S1_minus_S0__merged__L1")
    )
    
    for (pair in merged_pairs) {
      id_1 <- get_model_key("primary_family", "merged", family_name, pair[1])
      id_2 <- get_model_key("primary_family", "merged", family_name, pair[2])
      if (!is.null(cohort_draw_store[[id_1]]) && !is.null(cohort_draw_store[[id_2]])) {
        col_idx <- match(yr, cohort_draw_store[[id_1]]$years)
        contrast_rows[[length(contrast_rows) + 1L]] <- compute_selected_contrast(
          draw_mat_1 = cohort_draw_store[[id_1]]$mean_risk_draws[, col_idx, drop = FALSE],
          draw_mat_2 = cohort_draw_store[[id_2]]$mean_risk_draws[, col_idx, drop = FALSE],
          years = yr,
          direction_label = paste0(pair[3], "__", family_name)
        ) %>%
          dplyr::mutate(model_id_1 = id_1, model_id_2 = id_2)
      }
    }
    
    prior_pairs <- list(
      c("primary_family", "sa3_no_external", "main_minus_noexternal", "PNU", "B-L0"),
      c("sa1_robust", "sa3_no_external", "robust_minus_noexternal", "PNU", "B-L0"),
      c("primary_family", "sa3_no_external", "main_minus_noexternal", "SNU", "B-L0"),
      c("sa1_robust", "sa3_no_external", "robust_minus_noexternal", "SNU", "B-L0"),
      c("primary_family", "sa3_no_external", "main_minus_noexternal", "merged", "B-L0S0"),
      c("sa1_robust", "sa3_no_external", "robust_minus_noexternal", "merged", "B-L0S0"),
      c("primary_family", "sa3_no_external", "main_minus_noexternal", "merged", "B-L0S1"),
      c("sa1_robust", "sa3_no_external", "robust_minus_noexternal", "merged", "B-L0S1")
    )
    
    for (pair in prior_pairs) {
      id_1 <- get_model_key(pair[1], pair[4], family_name, pair[5])
      id_2 <- get_model_key(pair[2], pair[4], family_name, pair[5])
      if (!is.null(cohort_draw_store[[id_1]]) && !is.null(cohort_draw_store[[id_2]])) {
        col_idx <- match(yr, cohort_draw_store[[id_1]]$years)
        contrast_rows[[length(contrast_rows) + 1L]] <- compute_selected_contrast(
          draw_mat_1 = cohort_draw_store[[id_1]]$mean_risk_draws[, col_idx, drop = FALSE],
          draw_mat_2 = cohort_draw_store[[id_2]]$mean_risk_draws[, col_idx, drop = FALSE],
          years = yr,
          direction_label = paste0(pair[3], "__", pair[4], "__", pair[5], "__", family_name)
        ) %>%
          dplyr::mutate(model_id_1 = id_1, model_id_2 = id_2)
      }
    }
  }
}

bayes_posterior_contrasts <- bind_rows_or_template(
  contrast_rows,
  tibble::tibble(
    model_id_1 = character(),
    model_id_2 = character(),
    year = double(),
    contrast_label = character(),
    diff_mean = double(),
    diff_median = double(),
    diff_lcl = double(),
    diff_ucl = double(),
    p_diff_gt_0 = double()
  )
) %>%
  dplyr::relocate(model_id_1, model_id_2, .before = year)

# 🔴 Export: source-of-truth files and master RDS ===============================
analysis_spec <- tibble::tibble(
  data_path = data_path,
  export_path = export_path,
  active_bundle_ids = paste(active_bundle_ids, collapse = ";"),
  prediction_years = paste(prediction_years, collapse = ";"),
  selected_contrast_years = paste(selected_contrast_years, collapse = ";"),
  calibration_years = paste(calibration_years, collapse = ";"),
  chains_main = chains_main,
  warmup_main = warmup_main,
  iter_sampling_main = iter_sampling_main,
  adapt_delta_main = adapt_delta_main,
  max_treedepth_main = max_treedepth_main,
  chains_refit = chains_refit,
  warmup_refit = warmup_refit,
  iter_sampling_refit = iter_sampling_refit,
  adapt_delta_refit = adapt_delta_refit,
  max_treedepth_refit = max_treedepth_refit,
  prediction_draws_max = prediction_draws_max,
  tail_fit_absdiff_threshold = tail_fit_absdiff_threshold,
  prior_spec_rstan_version = "bayes_meta_informative_incidence_v_final_rstan_inline_v2",
  incidence_core_structure = "fixed informative incidence with sex + age_group (<20,20-29,30+)",
  latency_primary_structure = "continuous age_s + sex + optional age_s:sex and optional site(latency only)",
  family_set = "exponential;weibull;loglogistic;lognormal",
  model_time_unit = "years",
  raw_time_unit = "days",
  remission_main_rule = "treated_as_censoring",
  remission_ic_rule = "IPCW-weighted Bayesian cure model",
  remission_cr_rule = "Aalen-Johansen + Fine-Gray"
)

bayes_model_registry_csv <- model_registry_final %>%
  dplyr::left_join(bayes_weight_summary, by = "model_id")

exported_files <- c(
  safe_write_csv(analysis_spec, file.path(export_path, "analysis_spec_bayes.csv"), overwrite_existing_exports),
  safe_write_csv(branch_qc, file.path(export_path, "bayes_branch_qc.csv"), overwrite_existing_exports),
  safe_write_csv(km_yearly_summary, file.path(export_path, "km_yearly_summary.csv"), overwrite_existing_exports),
  safe_write_csv(bayes_model_registry_csv, file.path(export_path, "bayes_model_registry.csv"), overwrite_existing_exports),
  safe_write_csv(bayes_subject_static, file.path(export_path, "bayes_subject_static.csv"), overwrite_existing_exports),
  safe_write_csv(bayes_subject_predictions, file.path(export_path, "bayes_subject_predictions.csv"), overwrite_existing_exports),
  safe_write_csv(bayes_yearly_predictions, file.path(export_path, "bayes_yearly_predictions.csv"), overwrite_existing_exports),
  safe_write_csv(bayes_metrics, file.path(export_path, "bayes_metrics.csv"), overwrite_existing_exports),
  safe_write_csv(bayes_prior_utility, file.path(export_path, "bayes_prior_utility.csv"), overwrite_existing_exports),
  safe_write_csv(bayes_posterior_contrasts, file.path(export_path, "bayes_posterior_contrasts.csv"), overwrite_existing_exports),
  safe_write_csv(bayes_mcmc_diagnostics, file.path(export_path, "bayes_mcmc_diagnostics.csv"), overwrite_existing_exports),
  safe_write_csv(bayes_calibration, file.path(export_path, "bayes_calibration.csv"), overwrite_existing_exports),
  safe_write_csv(remission_cr_yearly, file.path(export_path, "bayes_remission_cr_yearly.csv"), overwrite_existing_exports),
  safe_write_csv(remission_crr_coef, file.path(export_path, "bayes_remission_crr_coef.csv"), overwrite_existing_exports)
)

master_results <- list(
  analysis_spec = analysis_spec,
  branch_qc = branch_qc,
  external_info = external_info,
  model_registry = model_registry_final,
  km_yearly_summary = km_yearly_summary,
  fits = fit_store,
  subject_static = bayes_subject_static,
  subject_predictions = bayes_subject_predictions,
  yearly_predictions = bayes_yearly_predictions,
  metrics = bayes_metrics,
  prior_utility = bayes_prior_utility,
  posterior_contrasts = bayes_posterior_contrasts,
  mcmc_diagnostics = bayes_mcmc_diagnostics,
  calibration = bayes_calibration,
  remission_competing_risk = list(
    yearly = remission_cr_yearly,
    finegray_coef = remission_crr_coef,
    raw_objects = remission_cr_results
  ),
  weight_summary = bayes_weight_summary,
  stan_code = stan_code
)

exported_files <- c(
  exported_files,
  safe_save_rds(master_results, file.path(export_path, "bayes_all_fits_master.rds"), overwrite_existing_exports)
)

invisible(exported_files)