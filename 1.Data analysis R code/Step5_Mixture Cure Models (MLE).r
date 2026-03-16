# 🔴 Configure: user-editable-paths-and-controls ===============================
data_path <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦5.Step5_MLE MCM/attachments'

## 🟠 Define: step5-analysis-controls ===============================
days_per_year <- 365.25
time_zero_epsilon <- 1e-06

run_datasets <- c("merged", "PNU", "SNU")
subgroup_var <- NULL

pred_times_years <- 1:10
cure_time_proxy_multiplier <- 50

event_code_transition <- 1L
censor_codes <- c(0L, 2L)

incidence_covariates <- c("site", "sex_fact", "age_exact_entry")
candidate_dists <- c("weibull", "lnorm", "llogis", "gengamma")
cure_link <- "logistic"
best_model_criterion <- "AIC"

# dist별 latency 공변량을 넣고 싶으면 아래 예시처럼 수정
# weibull / llogis: scale, shape
# lnorm: meanlog, sdlog
# gengamma: mu, sigma, Q
latency_anc_by_dist <- list(
  weibull = NULL,
  lnorm = NULL,
  llogis = NULL,
  gengamma = NULL
)
# 예시:
# latency_anc_by_dist$weibull <- list(scale = ~ age_exact_entry + sex_fact)
# latency_anc_by_dist$lnorm   <- list(meanlog = ~ age_exact_entry + sex_fact)
# latency_anc_by_dist$llogis  <- list(scale = ~ age_exact_entry + sex_fact)
# latency_anc_by_dist$gengamma <- list(mu = ~ age_exact_entry + sex_fact)

min_n_for_fit <- 30L
min_events_for_fit <- 10L
tail_start_year <- 5
tail_mae_good_cutoff <- 0.05

compute_time_dependent_auc <- TRUE
compute_ipcw_brier <- TRUE
compute_ibs <- TRUE
write_overlay_png <- TRUE

global_seed <- 20260301

# 🔴 Bootstrap: package-runtime-and-session-options ===============================
## 🟠 Ensure: required-packages-available ===============================
required_pkgs <- c(
  "survival", "flexsurv", "flexsurvcure", "dplyr", "tidyr", "purrr",
  "readr", "tibble", "stringr", "ggplot2", "timeROC"
)

ensure_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

dir.create(export_path, recursive = TRUE, showWarnings = FALSE)
set.seed(global_seed)
options(stringsAsFactors = FALSE)
ensure_packages(required_pkgs)

# 🔴 Declare: helper-functions-for-step5 ===============================
## 🟠 Build: data-qc-and-formula-helpers ===============================
safe_file_stub <- function(x) {
  x |>
    stringr::str_replace_all("[^[:alnum:]_\\-]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

stop_if_missing_columns <- function(data, cols) {
  missing_cols <- setdiff(cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("필수 컬럼이 없습니다: ", paste(missing_cols, collapse = ", "))
  }
}

drop_constant_terms <- function(data, terms_vec) {
  terms_vec <- terms_vec[terms_vec %in% names(data)]
  keep <- vapply(
    terms_vec,
    function(v) dplyr::n_distinct(stats::na.omit(data[[v]])) > 1,
    logical(1)
  )
  terms_vec[keep]
}

build_surv_formula <- function(rhs_terms) {
  rhs <- if (length(rhs_terms) == 0) "1" else paste(rhs_terms, collapse = " + ")
  stats::as.formula(
    paste0("survival::Surv(time_years_model, event_transition) ~ ", rhs)
  )
}

collect_model_vars <- function(formula_obj, anc_spec = NULL) {
  base_vars <- setdiff(all.vars(formula_obj), c("Surv"))
  anc_vars <- character(0)
  if (!is.null(anc_spec)) {
    anc_vars <- unique(unlist(lapply(anc_spec, all.vars), use.names = FALSE))
  }
  unique(c(base_vars, anc_vars))
}

prepare_analysis_data <- function(raw_df) {
  stop_if_missing_columns(raw_df, c("id", "site", "days_followup", "status_num"))
  
  dat <- raw_df |>
    dplyr::mutate(
      id = as.character(id),
      site = toupper(trimws(as.character(site))),
      sex_fact = dplyr::case_when(
        "sex_fact" %in% names(raw_df) ~ as.character(sex_fact),
        "sex_num" %in% names(raw_df) & suppressWarnings(as.numeric(sex_num)) == 0 ~ "Female",
        "sex_num" %in% names(raw_df) & suppressWarnings(as.numeric(sex_num)) == 1 ~ "Male",
        TRUE ~ NA_character_
      ),
      sex_fact = factor(sex_fact, levels = c("Female", "Male")),
      status_num = suppressWarnings(as.integer(status_num)),
      days_followup = suppressWarnings(as.numeric(days_followup)),
      age_exact_entry = if ("age_exact_entry" %in% names(raw_df)) suppressWarnings(as.numeric(age_exact_entry)) else NA_real_,
      age_int = if ("age_int" %in% names(raw_df)) suppressWarnings(as.numeric(age_int)) else NA_real_,
      event_transition = dplyr::if_else(status_num == event_code_transition, 1L, 0L, missing = 0L),
      time_years_raw = days_followup / days_per_year,
      time_years_model = pmax(time_years_raw, time_zero_epsilon),
      site_id = paste(site, id, sep = "::")
    ) |>
    dplyr::filter(!is.na(days_followup), !is.na(status_num), days_followup >= 0)
  
  if (!all(dat$status_num %in% c(0L, 1L, 2L))) {
    stop("status_num에는 0, 1, 2 외의 값이 포함되어 있습니다.")
  }
  
  if (nrow(dat) == 0) {
    stop("전처리 후 남은 데이터가 없습니다.")
  }
  
  dat
}

make_analysis_splits <- function(data, dataset_names, subgroup_var = NULL) {
  split_list <- list()
  idx <- 1L
  
  for (ds in dataset_names) {
    ds_upper <- toupper(ds)
    branch_df <- if (ds_upper == "MERGED") {
      data
    } else {
      dplyr::filter(data, site == ds_upper)
    }
    
    if (nrow(branch_df) == 0) next
    
    split_list[[idx]] <- list(
      dataset = ds,
      subgroup = "overall",
      data = branch_df
    )
    idx <- idx + 1L
    
    if (!is.null(subgroup_var)) {
      if (!subgroup_var %in% names(branch_df)) {
        stop("subgroup_var가 데이터에 없습니다: ", subgroup_var)
      }
      subgroup_vals <- unique(stats::na.omit(branch_df[[subgroup_var]]))
      subgroup_vals <- subgroup_vals[order(as.character(subgroup_vals))]
      
      for (sv in subgroup_vals) {
        sv_chr <- as.character(sv)
        sub_df <- branch_df |> dplyr::filter(.data[[subgroup_var]] == sv)
        if (nrow(sub_df) == 0) next
        
        split_list[[idx]] <- list(
          dataset = ds,
          subgroup = paste0(subgroup_var, "=", sv_chr),
          data = sub_df
        )
        idx <- idx + 1L
      }
    }
  }
  
  split_list
}

## 🟠 Compute: km-summary-and-ipcw-helpers ===============================
km_at_times <- function(time, event, times) {
  fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  sum_obj <- summary(fit, times = times, extend = TRUE)
  
  max_event_time <- if (any(event == 1L)) max(time[event == 1L]) else NA_real_
  max_followup_time <- max(time)
  
  tibble::tibble(
    year = times,
    km_surv = sum_obj$surv,
    km_risk = 1 - sum_obj$surv,
    km_lcl = sum_obj$lower,
    km_ucl = sum_obj$upper,
    n_risk = sum_obj$n.risk,
    max_event_time = max_event_time,
    max_followup_time = max_followup_time,
    support_flag = dplyr::case_when(
      is.na(max_event_time) ~ "no_events",
      times <= max_event_time ~ "observed_support",
      times <= max_followup_time ~ "beyond_last_event_censor_only",
      TRUE ~ "beyond_followup"
    )
  )
}

make_censor_survival_getters <- function(time, event) {
  censor_fit <- survival::survfit(survival::Surv(time, 1L - event) ~ 1)
  
  step_fn <- stats::stepfun(
    x = censor_fit$time,
    y = c(1, censor_fit$surv),
    right = TRUE
  )
  
  list(
    G_at = function(t) {
      pmax(step_fn(t), 1e-08)
    },
    G_left = function(t) {
      pmax(step_fn(pmax(t - 1e-10, 0)), 1e-08)
    }
  )
}

calc_ipcw_brier <- function(time, event, pred_surv, eval_time, G_at, G_left) {
  if (length(pred_surv) != length(time)) return(NA_real_)
  
  y_t <- as.integer(time > eval_time)
  w_t <- ifelse(
    time <= eval_time & event == 1L,
    1 / G_left(time),
    ifelse(time > eval_time, 1 / G_at(eval_time), 0)
  )
  
  mean(w_t * (y_t - pred_surv)^2, na.rm = TRUE)
}

calc_time_auc <- function(time, event, marker_risk, eval_time) {
  if (length(marker_risk) != length(time)) return(NA_real_)
  if (sum(event == 1L & time <= eval_time) == 0) return(NA_real_)
  if (sum(time > eval_time) == 0) return(NA_real_)
  if (length(unique(stats::na.omit(marker_risk))) < 2) return(NA_real_)
  
  auc_obj <- tryCatch(
    timeROC::timeROC(
      T = time,
      delta = event,
      marker = marker_risk,
      cause = 1,
      weighting = "marginal",
      times = eval_time,
      iid = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(auc_obj)) return(NA_real_)
  as.numeric(utils::tail(auc_obj$AUC, 1))
}

calc_trapezoid_ibs <- function(time_vec, brier_vec) {
  ok <- is.finite(time_vec) & is.finite(brier_vec)
  if (sum(ok) == 0) return(NA_real_)
  if (sum(ok) == 1) return(as.numeric(brier_vec[ok][1]))
  
  x <- time_vec[ok]
  y <- brier_vec[ok]
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  
  if (length(unique(x)) == 1) return(as.numeric(y[1]))
  
  area <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  area / (max(x) - min(x))
}

## 🟠 Predict: subject-level-and-mean-survival-helpers ===============================
predict_survival_long <- function(fit, newdata, times_vec) {
  pred_raw <- stats::predict(
    fit,
    newdata = newdata,
    type = "survival",
    times = times_vec,
    conf.int = FALSE
  ) |>
    tibble::as_tibble()
  
  pred_raw$.pred_row <- seq_len(nrow(pred_raw))
  
  if (".pred" %in% names(pred_raw)) {
    long_df <- pred_raw |>
      dplyr::select(.pred_row, .pred) |>
      tidyr::unnest(cols = .pred)
  } else {
    long_df <- pred_raw
  }
  
  time_col <- intersect(c(".time", ".eval_time", "time"), names(long_df))
  surv_col <- intersect(c(".pred_survival", "est", "survival"), names(long_df))
  
  if (length(time_col) == 0) {
    long_df[[".time"]] <- times_vec[1]
    time_col <- ".time"
  } else {
    time_col <- time_col[1]
  }
  
  if (length(surv_col) == 0) {
    stop("predict() 결과에서 survival probability 컬럼을 찾을 수 없습니다.")
  } else {
    surv_col <- surv_col[1]
  }
  
  long_df |>
    dplyr::transmute(
      .pred_row = .pred_row,
      year = as.numeric(.data[[time_col]]),
      surv_prob = as.numeric(.data[[surv_col]])
    )
}

## 🟠 Fit: one-distribution-mixture-cure-helper ===============================
fit_one_distribution <- function(sub_df, dataset_label, subgroup_label, dist_name, km_yearly_df) {
  fit_warnings <- character(0)
  
  incidence_terms_use <- drop_constant_terms(sub_df, incidence_covariates)
  fit_formula <- build_surv_formula(incidence_terms_use)
  anc_spec <- latency_anc_by_dist[[dist_name]]
  
  model_vars <- collect_model_vars(fit_formula, anc_spec)
  model_vars <- unique(c("site_id", "id", "site", model_vars))
  
  model_df <- sub_df |>
    dplyr::select(dplyr::any_of(unique(model_vars))) |>
    dplyr::mutate(.pred_row = dplyr::row_number()) |>
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::everything())))
  
  n_model <- nrow(model_df)
  n_events <- sum(model_df$event_transition, na.rm = TRUE)
  
  registry_stub <- tibble::tibble(
    dataset = dataset_label,
    subgroup = subgroup_label,
    model_class = "mixcure_MLE",
    dist = dist_name,
    n = n_model,
    events = n_events
  )
  
  if (n_model < min_n_for_fit) {
    return(list(
      registry = dplyr::mutate(
        registry_stub,
        fit_status = "skipped_low_n",
        logLik = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        criterion = NA_real_,
        warnings = NA_character_,
        error_message = paste0("n < ", min_n_for_fit)
      )
    ))
  }
  
  if (n_events < min_events_for_fit) {
    return(list(
      registry = dplyr::mutate(
        registry_stub,
        fit_status = "skipped_low_events",
        logLik = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        criterion = NA_real_,
        warnings = NA_character_,
        error_message = paste0("events < ", min_events_for_fit)
      )
    ))
  }
  
  fit_obj <- tryCatch(
    withCallingHandlers(
      {
        do.call(
          flexsurvcure::flexsurvcure,
          args = c(
            list(
              formula = fit_formula,
              data = model_df,
              dist = dist_name,
              link = cure_link,
              mixture = TRUE
            ),
            if (is.null(anc_spec)) list() else list(anc = anc_spec)
          )
        )
      },
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (inherits(fit_obj, "error")) {
    return(list(
      registry = dplyr::mutate(
        registry_stub,
        fit_status = "fit_failed",
        logLik = NA_real_,
        AIC = NA_real_,
        BIC = NA_real_,
        criterion = NA_real_,
        warnings = paste(unique(fit_warnings), collapse = " | "),
        error_message = conditionMessage(fit_obj)
      )
    ))
  }
  
  pred_long_core <- predict_survival_long(fit_obj, model_df, pred_times_years)
  
  pred_long <- pred_long_core |>
    dplyr::left_join(
      model_df |>
        dplyr::select(.pred_row, site_id, id, site),
      by = ".pred_row"
    ) |>
    dplyr::mutate(
      dataset = dataset_label,
      subgroup = subgroup_label,
      model_class = "mixcure_MLE",
      dist = dist_name,
      risk_prob = 1 - surv_prob,
      is_extrapolated = year > max(model_df$time_years_model),
      n = n_model,
      events = n_events
    ) |>
    dplyr::left_join(
      km_yearly_df |>
        dplyr::select(year, support_flag, km_surv, km_risk, n_risk),
      by = "year"
    )
  
  pred_mean <- pred_long |>
    dplyr::group_by(dataset, subgroup, model_class, dist, year) |>
    dplyr::summarise(
      mean_surv_prob = mean(surv_prob, na.rm = TRUE),
      mean_risk_prob = mean(risk_prob, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      km_yearly_df |>
        dplyr::select(year, km_surv, km_risk, support_flag, n_risk),
      by = "year"
    )
  
  censor_getters <- make_censor_survival_getters(
    time = model_df$time_years_model,
    event = model_df$event_transition
  )
  
  perf_by_time <- purrr::map_dfr(
    pred_times_years,
    function(tt) {
      pred_tt <- pred_long |>
        dplyr::filter(year == tt) |>
        dplyr::arrange(.pred_row)
      
      auc_val <- if (compute_time_dependent_auc) {
        calc_time_auc(
          time = model_df$time_years_model,
          event = model_df$event_transition,
          marker_risk = pred_tt$risk_prob,
          eval_time = tt
        )
      } else {
        NA_real_
      }
      
      brier_val <- if (compute_ipcw_brier) {
        calc_ipcw_brier(
          time = model_df$time_years_model,
          event = model_df$event_transition,
          pred_surv = pred_tt$surv_prob,
          eval_time = tt,
          G_at = censor_getters$G_at,
          G_left = censor_getters$G_left
        )
      } else {
        NA_real_
      }
      
      tibble::tibble(
        dataset = dataset_label,
        subgroup = subgroup_label,
        model_class = "mixcure_MLE",
        dist = dist_name,
        year = tt,
        AUC_t = auc_val,
        Brier_t = brier_val
      )
    }
  ) |>
    dplyr::left_join(
      km_yearly_df |>
        dplyr::select(year, support_flag, km_surv, km_risk, n_risk),
      by = "year"
    )
  
  ibs_val <- if (compute_ibs && compute_ipcw_brier) {
    calc_trapezoid_ibs(perf_by_time$year, perf_by_time$Brier_t)
  } else {
    NA_real_
  }
  
  cure_proxy_time <- max(
    pred_times_years,
    model_df$time_years_model,
    na.rm = TRUE
  ) * cure_time_proxy_multiplier
  
  cure_detail <- predict_survival_long(fit_obj, model_df, cure_proxy_time) |>
    dplyr::left_join(
      model_df |>
        dplyr::select(.pred_row, site_id, id, site),
      by = ".pred_row"
    ) |>
    dplyr::transmute(
      dataset = dataset_label,
      subgroup = subgroup_label,
      model_class = "mixcure_MLE",
      dist = dist_name,
      site_id = site_id,
      id = id,
      site = site,
      cure_proxy_time_years = year,
      approx_cure_fraction = surv_prob
    )
  
  tail_df <- pred_mean |>
    dplyr::filter(
      year >= tail_start_year,
      support_flag != "beyond_followup"
    )
  
  tail_mae <- if (nrow(tail_df) > 0) {
    mean(abs(tail_df$mean_surv_prob - tail_df$km_surv), na.rm = TRUE)
  } else {
    NA_real_
  }
  
  tail_fit_flag <- dplyr::case_when(
    is.na(tail_mae) ~ NA_character_,
    tail_mae <= tail_mae_good_cutoff ~ "good",
    tail_mae <= 2 * tail_mae_good_cutoff ~ "moderate",
    TRUE ~ "poor"
  )
  
  coef_tbl <- tibble::tibble(
    dataset = dataset_label,
    subgroup = subgroup_label,
    model_class = "mixcure_MLE",
    dist = dist_name,
    term = names(stats::coef(fit_obj)),
    estimate = unname(stats::coef(fit_obj))
  )
  
  perf_model_level <- tibble::tibble(
    dataset = dataset_label,
    subgroup = subgroup_label,
    model_class = "mixcure_MLE",
    dist = dist_name,
    n = n_model,
    events = n_events,
    logLik = as.numeric(stats::logLik(fit_obj)),
    AIC = stats::AIC(fit_obj),
    BIC = stats::BIC(fit_obj),
    IBS = ibs_val,
    approx_cure_fraction_mean = mean(cure_detail$approx_cure_fraction, na.rm = TRUE),
    tail_mae = tail_mae,
    tail_fit_flag = tail_fit_flag,
    fit_status = "success",
    warnings = paste(unique(fit_warnings), collapse = " | "),
    error_message = NA_character_
  ) |>
    dplyr::mutate(
      criterion = dplyr::case_when(
        best_model_criterion == "AIC" ~ AIC,
        best_model_criterion == "BIC" ~ BIC,
        TRUE ~ AIC
      )
    )
  
  fit_key <- safe_file_stub(paste(dataset_label, subgroup_label, dist_name, sep = "__"))
  
  list(
    registry = perf_model_level,
    fit_key = fit_key,
    fit = fit_obj,
    pred_subject = pred_long,
    pred_mean = pred_mean,
    perf_by_time = perf_by_time,
    cure_detail = cure_detail,
    coef_tbl = coef_tbl
  )
}

# 🔴 Read: source-data-and-build-step5-input ===============================
## 🟠 Import: merged-dataset3-from-csv ===============================
raw_data <- readr::read_csv(
  file = data_path,
  show_col_types = FALSE,
  guess_max = 100000
)

## 🟠 Transform: transition-focused-analysis-data ===============================
dat_analysis <- prepare_analysis_data(raw_data)

## 🟠 Check: minimal-event-and-site-coverage ===============================
if (!any(dat_analysis$event_transition == 1L)) {
  stop("transition event가 하나도 없습니다. Step5를 진행할 수 없습니다.")
}

if (!all(toupper(run_datasets) %in% c("MERGED", unique(dat_analysis$site)))) {
  warning("run_datasets 중 일부가 데이터에 존재하지 않을 수 있습니다: ",
          paste(setdiff(toupper(run_datasets), c("MERGED", unique(dat_analysis$site))), collapse = ", "))
}

saveRDS(
  object = dat_analysis,
  file = file.path(export_path, "step5_analysis_input.rds")
)

# 🔴 Construct: dataset-and-subgroup-analysis-plan ===============================
## 🟠 Create: branch-specific-split-list ===============================
analysis_splits <- make_analysis_splits(
  data = dat_analysis,
  dataset_names = run_datasets,
  subgroup_var = subgroup_var
)

if (length(analysis_splits) == 0) {
  stop("분석 가능한 dataset/subgroup split이 없습니다.")
}

km_yearly_all <- purrr::map_dfr(
  analysis_splits,
  function(x) {
    km_at_times(
      time = x$data$time_years_model,
      event = x$data$event_transition,
      times = pred_times_years
    ) |>
      dplyr::mutate(
        dataset = x$dataset,
        subgroup = x$subgroup,
        model_class = "plain_KM",
        dist = NA_character_
      ) |>
      dplyr::select(
        dataset, subgroup, model_class, dist, year,
        km_surv, km_risk, km_lcl, km_ucl, n_risk,
        max_event_time, max_followup_time, support_flag
      )
  }
)

# 🔴 Execute: mixture-cure-mle-fitting-across-splits ===============================
## 🟠 Iterate: dataset-subgroup-and-latency-family ===============================
fit_results_nested <- list()
fit_objects <- list()
pred_subject_long <- list()
pred_mean_by_time <- list()
perf_by_time_all <- list()
perf_model_level_all <- list()
cure_detail_all <- list()
coef_table_all <- list()

result_idx <- 1L

for (split_obj in analysis_splits) {
  ds_label <- split_obj$dataset
  sg_label <- split_obj$subgroup
  sub_df <- split_obj$data
  
  km_one <- km_yearly_all |>
    dplyr::filter(dataset == ds_label, subgroup == sg_label)
  
  for (dist_name in candidate_dists) {
    message("Fitting Step5 mixture cure model: dataset=", ds_label,
            ", subgroup=", sg_label, ", dist=", dist_name)
    
    one_res <- fit_one_distribution(
      sub_df = sub_df,
      dataset_label = ds_label,
      subgroup_label = sg_label,
      dist_name = dist_name,
      km_yearly_df = km_one
    )
    
    fit_results_nested[[result_idx]] <- one_res
    result_idx <- result_idx + 1L
    
    perf_model_level_all[[length(perf_model_level_all) + 1L]] <- one_res$registry
    
    if (!is.null(one_res$fit)) {
      fit_objects[[one_res$fit_key]] <- one_res$fit
      pred_subject_long[[length(pred_subject_long) + 1L]] <- one_res$pred_subject
      pred_mean_by_time[[length(pred_mean_by_time) + 1L]] <- one_res$pred_mean
      perf_by_time_all[[length(perf_by_time_all) + 1L]] <- one_res$perf_by_time
      cure_detail_all[[length(cure_detail_all) + 1L]] <- one_res$cure_detail
      coef_table_all[[length(coef_table_all) + 1L]] <- one_res$coef_tbl
    }
  }
}

fit_registry <- dplyr::bind_rows(perf_model_level_all) |>
  dplyr::mutate(
    fit_key = dplyr::if_else(
      fit_status == "success",
      safe_file_stub(paste(dataset, subgroup, dist, sep = "__")),
      NA_character_
    )
  ) |>
  dplyr::arrange(dataset, subgroup, dplyr::across(dplyr::all_of(best_model_criterion)))

pred_subject_long_df <- dplyr::bind_rows(pred_subject_long)
pred_mean_by_time_df <- dplyr::bind_rows(pred_mean_by_time)
perf_by_time_df <- dplyr::bind_rows(perf_by_time_all)
cure_detail_df <- dplyr::bind_rows(cure_detail_all)
coef_table_df <- dplyr::bind_rows(coef_table_all)

# 🔴 Select: best-fitting-step5-model-per-analysis-unit ===============================
## 🟠 Rank: candidate-models-by-chosen-criterion ===============================
best_model_registry <- fit_registry |>
  dplyr::filter(fit_status == "success") |>
  dplyr::group_by(dataset, subgroup) |>
  dplyr::slice_min(order_by = .data[[best_model_criterion]], n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

best_fit_keys <- stats::setNames(best_model_registry$fit_key, paste(best_model_registry$dataset, best_model_registry$subgroup, sep = "__"))
best_fit_objects <- fit_objects[best_model_registry$fit_key]

best_pred_mean_by_time_df <- pred_mean_by_time_df |>
  dplyr::inner_join(
    best_model_registry |>
      dplyr::select(dataset, subgroup, dist),
    by = c("dataset", "subgroup", "dist")
  )

best_pred_subject_long_df <- pred_subject_long_df |>
  dplyr::inner_join(
    best_model_registry |>
      dplyr::select(dataset, subgroup, dist),
    by = c("dataset", "subgroup", "dist")
  )

best_perf_by_time_df <- perf_by_time_df |>
  dplyr::inner_join(
    best_model_registry |>
      dplyr::select(dataset, subgroup, dist),
    by = c("dataset", "subgroup", "dist")
  )

best_cure_detail_df <- cure_detail_df |>
  dplyr::inner_join(
    best_model_registry |>
      dplyr::select(dataset, subgroup, dist),
    by = c("dataset", "subgroup", "dist")
  )

best_coef_table_df <- coef_table_df |>
  dplyr::inner_join(
    best_model_registry |>
      dplyr::select(dataset, subgroup, dist),
    by = c("dataset", "subgroup", "dist")
  )

# 🔴 Assemble: overlay-data-from-source-of-truth-dataframes ===============================
## 🟠 Merge: km-and-model-curves-into-plotting-dataframe ===============================
overlay_curve_data <- dplyr::bind_rows(
  km_yearly_all |>
    dplyr::transmute(
      dataset = dataset,
      subgroup = subgroup,
      curve_source = "plain_KM",
      model_class = "plain_KM",
      dist = NA_character_,
      year = year,
      surv_value = km_surv,
      risk_value = km_risk,
      lower_value = km_lcl,
      upper_value = km_ucl,
      support_flag = support_flag
    ),
  pred_mean_by_time_df |>
    dplyr::transmute(
      dataset = dataset,
      subgroup = subgroup,
      curve_source = "mixcure_MLE",
      model_class = model_class,
      dist = dist,
      year = year,
      surv_value = mean_surv_prob,
      risk_value = mean_risk_prob,
      lower_value = NA_real_,
      upper_value = NA_real_,
      support_flag = support_flag
    )
)

best_overlay_curve_data <- dplyr::bind_rows(
  km_yearly_all |>
    dplyr::transmute(
      dataset = dataset,
      subgroup = subgroup,
      curve_source = "plain_KM",
      model_class = "plain_KM",
      dist = NA_character_,
      year = year,
      surv_value = km_surv,
      risk_value = km_risk,
      lower_value = km_lcl,
      upper_value = km_ucl,
      support_flag = support_flag
    ),
  best_pred_mean_by_time_df |>
    dplyr::transmute(
      dataset = dataset,
      subgroup = subgroup,
      curve_source = paste0("mixcure_MLE_", dist),
      model_class = model_class,
      dist = dist,
      year = year,
      surv_value = mean_surv_prob,
      risk_value = mean_risk_prob,
      lower_value = NA_real_,
      upper_value = NA_real_,
      support_flag = support_flag
    )
)

# 🔴 Save: step5-rds-bundles-for-reuse ===============================
## 🟠 Persist: all-candidate-and-best-fit-objects ===============================
step5_fit_bundle_all <- list(
  config = list(
    data_path = data_path,
    export_path = export_path,
    run_datasets = run_datasets,
    subgroup_var = subgroup_var,
    pred_times_years = pred_times_years,
    incidence_covariates = incidence_covariates,
    candidate_dists = candidate_dists,
    cure_link = cure_link,
    best_model_criterion = best_model_criterion,
    latency_anc_by_dist = latency_anc_by_dist
  ),
  analysis_data = dat_analysis,
  fit_registry = fit_registry,
  fit_objects = fit_objects,
  km_yearly_all = km_yearly_all,
  pred_mean_by_time = pred_mean_by_time_df,
  perf_by_time = perf_by_time_df,
  perf_model_level = fit_registry,
  cure_detail = cure_detail_df,
  coef_table = coef_table_df
)

step5_fit_bundle_best <- list(
  config = step5_fit_bundle_all$config,
  analysis_data = dat_analysis,
  best_model_registry = best_model_registry,
  best_fit_objects = best_fit_objects,
  km_yearly_all = km_yearly_all,
  pred_subject_long = best_pred_subject_long_df,
  pred_mean_by_time = best_pred_mean_by_time_df,
  perf_by_time = best_perf_by_time_df,
  cure_detail = best_cure_detail_df,
  coef_table = best_coef_table_df
)

saveRDS(
  object = step5_fit_bundle_all,
  file = file.path(export_path, "step5_fit_bundle_all.rds")
)

saveRDS(
  object = step5_fit_bundle_best,
  file = file.path(export_path, "step5_fit_bundle_best.rds")
)

# 🔴 Write: csv-source-of-truth-exports ===============================
## 🟠 Export: model-registry-and-performance-csv ===============================
readr::write_csv(
  fit_registry,
  file.path(export_path, "step5_fit_registry.csv")
)

readr::write_csv(
  pred_mean_by_time_df,
  file.path(export_path, "step5_pred_mean_by_time.csv")
)

readr::write_csv(
  perf_by_time_df,
  file.path(export_path, "step5_perf_by_time.csv")
)

readr::write_csv(
  fit_registry |>
    dplyr::select(
      dataset, subgroup, model_class, dist, n, events,
      logLik, AIC, BIC, IBS, approx_cure_fraction_mean,
      tail_mae, tail_fit_flag, fit_status, warnings, error_message, criterion, fit_key
    ),
  file.path(export_path, "step5_perf_model_level.csv")
)

## 🟠 Export: subject-level-and-cure-detail-csv ===============================
readr::write_csv(
  pred_subject_long_df |>
    dplyr::select(
      dataset, subgroup, model_class, dist, site_id, id, site,
      year, surv_prob, risk_prob, is_extrapolated,
      support_flag, km_surv, km_risk, n_risk
    ),
  file.path(export_path, "step5_pred_subject_long.csv.gz")
)

readr::write_csv(
  cure_detail_df,
  file.path(export_path, "step5_cure_detail_subject.csv.gz")
)

readr::write_csv(
  coef_table_df,
  file.path(export_path, "step5_coef_table.csv")
)

## 🟠 Export: km-overlay-dataframe-csv ===============================
readr::write_csv(
  overlay_curve_data,
  file.path(export_path, "step5_overlay_curve_data.csv")
)

# 🔴 Render: overlay-png-from-exported-dataframe ===============================
## 🟠 Draw: best-model-versus-km-faceted-overlay ===============================
if (write_overlay_png && nrow(best_overlay_curve_data) > 0) {
  overlay_plot <- ggplot2::ggplot(
    best_overlay_curve_data,
    ggplot2::aes(
      x = year,
      y = surv_value,
      color = curve_source,
      linetype = curve_source
    )
  ) +
    ggplot2::geom_step(linewidth = 0.7) +
    ggplot2::facet_grid(subgroup ~ dataset, scales = "fixed") +
    ggplot2::scale_x_continuous(breaks = pred_times_years) +
    ggplot2::labs(
      x = "Years since cohort entry",
      y = "Overall survival probability",
      color = "Curve",
      linetype = "Curve",
      title = "Step5: Best frequentist mixture cure model vs plain KM"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  ggplot2::ggsave(
    filename = file.path(export_path, "step5_overlay_best_models.png"),
    plot = overlay_plot,
    width = 12,
    height = max(4, 2 + 2 * max(1, dplyr::n_distinct(best_overlay_curve_data$subgroup))),
    dpi = 300
  )
}