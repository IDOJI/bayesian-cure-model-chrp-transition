# 🔴 Configure: paths, horizons, and file budget ===============================

## 🟠 Declare: user directories and reporting times ===============================
path_combined <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧0.생존 데이터 처리와 요약/🟦2.데이터3 처리/attachments/MERGED_dataset3_pnu_snu.csv"
export_dir <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦1.Step1_KM/attachments"

eval_years <- 1:10
eval_times_days <- eval_years * 365.25

strict_schema <- TRUE
plot_width_in <- 8
plot_height_in <- 6
plot_dpi <- 300
max_generated_files <- 20L

STATUS_LEVELS <- c("right_censoring", "remission", "transition")
EXPECTED_SITES <- c("PNU", "SNU")

## 🟠 Attach: package set and export cap ===============================
required_pkgs <- c("survival", "readr", "dplyr", "tibble", "ggplot2")

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# 🔴 Define: schema guards and export constructors ===============================

## 🟠 Implement: scalar helpers and label mapping ===============================
write_csv_safe <- function(x, path) {
  readr::write_csv(x, path, na = "")
}

fill_na_with <- function(x, fill) {
  out <- x
  out[is.na(out)] <- fill[is.na(out)]
  out
}

safe_num <- function(x, name) {
  if (name %in% names(x)) as.numeric(unname(x[[name]])) else NA_real_
}

as_scalar_character <- function(x, name) {
  x <- as.character(x)
  if (length(x) < 1L || is.na(x[1]) || x[1] == "") {
    stop(name, " must be a non-missing scalar character value.")
  }
  x[1]
}

as_scalar_numeric <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) < 1L || !is.finite(x[1])) {
    stop(name, " must be a finite scalar numeric value.")
  }
  x[1]
}

status_num_to_label <- function(x) {
  dplyr::case_when(
    x == 0L ~ "right_censoring",
    x == 1L ~ "transition",
    x == 2L ~ "remission",
    TRUE    ~ NA_character_
  )
}

standardize_status_factor <- function(status_num) {
  factor(status_num_to_label(status_num), levels = STATUS_LEVELS)
}

sanitize_prefix <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

as_date_if_present <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.character(x)) return(as.Date(x))
  x
}

## 🟠 Implement: merged-data validation rules ===============================
run_basic_qc <- function(dat, strict_schema = TRUE) {
  required_cols <- c("id", "site", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  dat <- dat |>
    dplyr::mutate(
      id = trimws(as.character(id)),
      site = toupper(trimws(as.character(site))),
      days_followup = as.numeric(days_followup),
      status_num = as.integer(status_num)
    )
  
  if (any(is.na(dat$id) | dat$id == "")) {
    stop("id contains missing or blank values.")
  }
  if (any(is.na(dat$site) | dat$site == "")) {
    stop("site contains missing or blank values.")
  }
  if (any(is.na(dat$days_followup))) {
    stop("days_followup contains NA.")
  }
  if (any(dat$days_followup < 0)) {
    stop("days_followup contains values < 0.")
  }
  if (any(is.na(dat$status_num))) {
    stop("status_num contains NA.")
  }
  if (!all(dat$status_num %in% c(0L, 1L, 2L))) {
    stop("status_num must be one of 0, 1, 2.")
  }
  
  dat <- dat |>
    dplyr::mutate(
      site_id = paste(site, id, sep = "_")
    )
  
  if (anyDuplicated(dat$site_id) > 0) {
    dup_keys <- unique(dat$site_id[duplicated(dat$site_id)])
    stop(
      "Duplicated site_id detected. Step 1 assumes one row per subject. Examples: ",
      paste(head(dup_keys, 5), collapse = ", ")
    )
  }
  
  expected_status <- status_num_to_label(dat$status_num)
  if ("status" %in% names(dat)) {
    observed_status <- trimws(as.character(dat$status))
    observed_status[observed_status == ""] <- NA_character_
    
    mismatch_n <- sum(
      !is.na(expected_status) &
        !is.na(observed_status) &
        expected_status != observed_status
    )
    
    if (mismatch_n > 0) {
      warning(
        "status and status_num are inconsistent in ", mismatch_n,
        " rows. status will be overwritten using status_num."
      )
    }
  }
  
  dat <- dat |>
    dplyr::mutate(
      status = standardize_status_factor(status_num)
    )
  
  if ("sex_num" %in% names(dat)) {
    dat <- dat |>
      dplyr::mutate(
        sex_num = dplyr::if_else(
          sex_num %in% c(0, 1),
          as.integer(sex_num),
          NA_integer_
        )
      )
  }
  
  if ("sex_fact" %in% names(dat)) {
    dat <- dat |>
      dplyr::mutate(
        sex_fact = factor(as.character(sex_fact), levels = c("Female", "Male"))
      )
  } else if ("sex_num" %in% names(dat)) {
    dat <- dat |>
      dplyr::mutate(
        sex_fact = dplyr::case_when(
          sex_num == 0L ~ "Female",
          sex_num == 1L ~ "Male",
          TRUE ~ NA_character_
        ),
        sex_fact = factor(sex_fact, levels = c("Female", "Male"))
      )
  }
  
  for (nm in intersect(c("date_birth", "date_entry"), names(dat))) {
    dat[[nm]] <- as_date_if_present(dat[[nm]])
  }
  
  unexpected_sites <- setdiff(sort(unique(dat$site)), EXPECTED_SITES)
  if (length(unexpected_sites) > 0) {
    warning(
      "Unexpected site values detected: ",
      paste(unexpected_sites, collapse = ", "),
      ". The script will still run subject to the file budget."
    )
  }
  
  snu_remission_n <- sum(dat$site == "SNU" & dat$status_num == 2L, na.rm = TRUE)
  if (snu_remission_n > 0) {
    msg <- paste0(
      "SNU has remission events (status_num == 2) in ", snu_remission_n,
      " rows. This violates the declared data dictionary."
    )
    if (isTRUE(strict_schema)) {
      stop(msg)
    } else {
      warning(msg)
    }
  }
  
  dat
}

## 🟠 Implement: dataset splitting and file targets ===============================
make_dataset_specs <- function(dat) {
  site_split <- split(dat, dat$site)
  site_split <- site_split[order(names(site_split))]
  
  specs <- list(
    list(
      dataset_label = "MERGED",
      file_prefix = "merged",
      data = dat
    )
  )
  
  for (site_name in names(site_split)) {
    specs[[length(specs) + 1L]] <- list(
      dataset_label = site_name,
      file_prefix = sanitize_prefix(site_name),
      data = site_split[[site_name]]
    )
  }
  
  prefixes <- vapply(specs, function(x) x$file_prefix, character(1))
  if (anyDuplicated(prefixes) > 0) {
    stop("Sanitized file_prefix values are duplicated. Please fix site labels.")
  }
  
  names(specs) <- prefixes
  specs
}

make_target_files <- function(dataset_specs, export_dir) {
  csv_files <- c(
    "step1_dataset_summary.csv",
    "step1_curve_source.csv",
    "step1_km_timepoints_primary.csv",
    "step1_km_timepoints_all.csv"
  )
  
  prefixes <- vapply(dataset_specs, function(x) x$file_prefix, character(1))
  
  png_files <- unlist(lapply(prefixes, function(px) {
    c(
      paste0(px, "_step1_survival_transitionfree.png"),
      paste0(px, "_step1_one_minus_km_transition.png"),
      paste0(px, "_step1_reversekm_followup.png")
    )
  }), use.names = FALSE)
  
  rds_files <- paste0(prefixes, "_step1_fit_bundle.rds")
  
  all_files <- c(csv_files, png_files, rds_files)
  file.path(export_dir, all_files)
}

## 🟠 Implement: survival table builders ===============================
build_curve_table <- function(fit, dataset_label, file_prefix, curve_type, max_followup_days) {
  dataset_label_scalar <- as_scalar_character(dataset_label, "dataset_label")
  file_prefix_scalar <- as_scalar_character(file_prefix, "file_prefix")
  curve_type_scalar <- as_scalar_character(curve_type, "curve_type")
  max_fu_days_scalar <- as_scalar_numeric(max_followup_days, "max_followup_days")
  max_fu_years_scalar <- max_fu_days_scalar / 365.25
  
  s <- summary(fit, censored = TRUE)
  
  if (length(s$time) == 0L) {
    n_rows <- 1L
    time_days_vec <- 0
    time_years_vec <- 0
    n_risk_vec <- as.numeric(fit$n)[1]
    n_event_vec <- 0
    n_censor_vec <- 0
    surv_vec <- 1
    se_vec <- 0
    lower_vec <- 1
    upper_vec <- 1
  } else {
    time_days_vec <- c(0, as.numeric(s$time))
    time_years_vec <- time_days_vec / 365.25
    surv_vec <- c(1, as.numeric(s$surv))
    se_vec <- c(0, fill_na_with(as.numeric(s$std.err), rep(0, length(s$surv))))
    lower_vec <- c(1, fill_na_with(as.numeric(s$lower), as.numeric(s$surv)))
    upper_vec <- c(1, fill_na_with(as.numeric(s$upper), as.numeric(s$surv)))
    n_risk_vec <- c(as.numeric(fit$n)[1], as.numeric(s$n.risk))
    n_event_vec <- c(0, as.numeric(s$n.event))
    n_censor_vec <- c(0, as.numeric(s$n.censor))
    n_rows <- length(time_days_vec)
    
    lengths_to_check <- c(
      length(time_days_vec), length(time_years_vec), length(surv_vec),
      length(se_vec), length(lower_vec), length(upper_vec),
      length(n_risk_vec), length(n_event_vec), length(n_censor_vec)
    )
    
    if (length(unique(lengths_to_check)) != 1L) {
      stop("Inconsistent summary vector lengths detected in build_curve_table().")
    }
  }
  
  out <- tibble::tibble(
    dataset_label = rep(dataset_label_scalar, n_rows),
    file_prefix = rep(file_prefix_scalar, n_rows),
    curve_type = rep(curve_type_scalar, n_rows),
    time_days = time_days_vec,
    time_years = time_years_vec,
    n_risk = n_risk_vec,
    n_event = n_event_vec,
    n_censor = n_censor_vec,
    surv = surv_vec,
    std_err = se_vec,
    lower_95 = lower_vec,
    upper_95 = upper_vec,
    max_followup_days = rep(max_fu_days_scalar, n_rows),
    max_followup_years = rep(max_fu_years_scalar, n_rows),
    within_observed_followup = time_days_vec <= max_fu_days_scalar
  )
  
  if (identical(curve_type_scalar, "km_transition")) {
    out <- out |>
      dplyr::mutate(
        one_minus_km_transition = 1 - surv,
        one_minus_km_transition_lcl95 = 1 - upper_95,
        one_minus_km_transition_ucl95 = 1 - lower_95
      )
  } else {
    out <- out |>
      dplyr::mutate(
        one_minus_km_transition = NA_real_,
        one_minus_km_transition_lcl95 = NA_real_,
        one_minus_km_transition_ucl95 = NA_real_
      )
  }
  
  out
}

build_timepoint_table <- function(fit, times_days, dataset_label, file_prefix, max_followup_days) {
  dataset_label_scalar <- as_scalar_character(dataset_label, "dataset_label")
  file_prefix_scalar <- as_scalar_character(file_prefix, "file_prefix")
  max_fu_days_scalar <- as_scalar_numeric(max_followup_days, "max_followup_days")
  max_fu_years_scalar <- max_fu_days_scalar / 365.25
  times_days <- as.numeric(times_days)
  
  s <- summary(fit, times = times_days, extend = TRUE)
  
  n_rows <- length(times_days)
  surv_vec <- as.numeric(s$surv)
  se_vec <- fill_na_with(as.numeric(s$std.err), rep(0, n_rows))
  lower_vec <- fill_na_with(as.numeric(s$lower), surv_vec)
  upper_vec <- fill_na_with(as.numeric(s$upper), surv_vec)
  n_risk_vec <- as.numeric(s$n.risk)
  n_event_vec <- as.numeric(s$n.event)
  n_censor_vec <- as.numeric(s$n.censor)
  
  lengths_to_check <- c(
    length(times_days), length(surv_vec), length(se_vec), length(lower_vec),
    length(upper_vec), length(n_risk_vec), length(n_event_vec), length(n_censor_vec)
  )
  
  if (length(unique(lengths_to_check)) != 1L) {
    stop("Inconsistent summary vector lengths detected in build_timepoint_table().")
  }
  
  tibble::tibble(
    dataset_label = rep(dataset_label_scalar, n_rows),
    file_prefix = rep(file_prefix_scalar, n_rows),
    curve_type = rep("km_transition", n_rows),
    time_days = times_days,
    time_years = times_days / 365.25,
    n_risk = n_risk_vec,
    n_event = n_event_vec,
    n_censor = n_censor_vec,
    surv = surv_vec,
    std_err = se_vec,
    lower_95 = lower_vec,
    upper_95 = upper_vec,
    one_minus_km_transition = 1 - surv_vec,
    one_minus_km_transition_lcl95 = 1 - upper_vec,
    one_minus_km_transition_ucl95 = 1 - lower_vec,
    max_followup_days = rep(max_fu_days_scalar, n_rows),
    max_followup_years = rep(max_fu_years_scalar, n_rows),
    within_observed_followup = times_days <= max_fu_days_scalar
  )
}

build_summary_row <- function(dat, km_fit, reverse_km, dataset_label, file_prefix) {
  dataset_label_scalar <- as_scalar_character(dataset_label, "dataset_label")
  file_prefix_scalar <- as_scalar_character(file_prefix, "file_prefix")
  
  dat2 <- dat |>
    dplyr::mutate(event_transition = as.integer(status_num == 1L))
  
  km_s <- summary(km_fit)
  rkm_tbl <- summary(reverse_km)$table
  max_fu <- max(dat2$days_followup, na.rm = TRUE)
  
  if (any(dat2$event_transition == 1L, na.rm = TRUE)) {
    last_transition <- max(dat2$days_followup[dat2$event_transition == 1L], na.rm = TRUE)
    idx <- max(which(km_s$time <= last_transition))
    n_risk_last <- as.numeric(km_s$n.risk[idx])
    surv_last <- as.numeric(km_s$surv[idx])
    
    plateau_length_days <- max_fu - last_transition
    one_more_event_drop_pct <- ifelse(
      is.na(n_risk_last) || n_risk_last <= 0,
      NA_real_,
      100 * surv_last / n_risk_last
    )
    n_right_censor_after_last_transition <- sum(
      dat2$status_num == 0L & dat2$days_followup > last_transition,
      na.rm = TRUE
    )
    n_remission_after_last_transition <- sum(
      dat2$status_num == 2L & dat2$days_followup > last_transition,
      na.rm = TRUE
    )
    n_primary_censor_after_last_transition <- sum(
      dat2$status_num != 1L & dat2$days_followup > last_transition,
      na.rm = TRUE
    )
  } else {
    last_transition <- NA_real_
    plateau_length_days <- NA_real_
    n_risk_last <- NA_real_
    surv_last <- NA_real_
    one_more_event_drop_pct <- NA_real_
    n_right_censor_after_last_transition <- NA_real_
    n_remission_after_last_transition <- NA_real_
    n_primary_censor_after_last_transition <- NA_real_
  }
  
  tibble::tibble(
    dataset_label = dataset_label_scalar,
    file_prefix = file_prefix_scalar,
    n = nrow(dat2),
    n_right_censoring = sum(dat2$status_num == 0L, na.rm = TRUE),
    n_transition = sum(dat2$status_num == 1L, na.rm = TRUE),
    n_remission = sum(dat2$status_num == 2L, na.rm = TRUE),
    n_primary_events = sum(dat2$status_num == 1L, na.rm = TRUE),
    n_primary_censored = sum(dat2$status_num != 1L, na.rm = TRUE),
    min_followup_days = min(dat2$days_followup, na.rm = TRUE),
    median_observed_time_days_raw = stats::median(dat2$days_followup, na.rm = TRUE),
    max_followup_days = max_fu,
    max_followup_years = max_fu / 365.25,
    reversekm_median_days = safe_num(rkm_tbl, "median"),
    reversekm_median_lcl95_days = safe_num(rkm_tbl, "0.95LCL"),
    reversekm_median_ucl95_days = safe_num(rkm_tbl, "0.95UCL"),
    reversekm_rmean_days = safe_num(rkm_tbl, "rmean"),
    reversekm_se_rmean_days = safe_num(rkm_tbl, "se(rmean)"),
    last_transition_days = last_transition,
    plateau_length_days = plateau_length_days,
    n_risk_at_last_transition = n_risk_last,
    km_at_last_transition = surv_last,
    one_more_event_drop_pct = one_more_event_drop_pct,
    n_right_censor_after_last_transition = n_right_censor_after_last_transition,
    n_remission_after_last_transition = n_remission_after_last_transition,
    n_primary_censor_after_last_transition = n_primary_censor_after_last_transition,
    one_minus_km_interpretation = "Descriptive 1-KM under remission-as-censoring assumption; NOT a CIF",
    reverse_km_interpretation = "Descriptive only; NOT a formal sufficient-follow-up test",
    reporting_rule = "Use step1_km_timepoints_primary.csv for main reporting; step1_km_timepoints_all.csv is supplementary only",
    betensky_stability_limits_implemented = FALSE
  )
}

## 🟠 Implement: per-dataset fitting bundle ===============================
fit_step1_dataset <- function(dat, dataset_label, file_prefix, eval_years, eval_times_days) {
  if (nrow(dat) == 0L) {
    stop("Dataset is empty: ", dataset_label)
  }
  
  dat <- dat |>
    dplyr::mutate(
      event_transition = as.integer(status_num == 1L),
      censoring_primary = as.integer(status_num != 1L)
    ) |>
    dplyr::arrange(days_followup)
  
  max_fu <- max(dat$days_followup, na.rm = TRUE)
  
  km_transition <- survival::survfit(
    survival::Surv(days_followup, event_transition) ~ 1,
    data = dat,
    conf.type = "log-log"
  )
  
  reverse_km <- survival::survfit(
    survival::Surv(days_followup, censoring_primary) ~ 1,
    data = dat,
    conf.type = "log-log"
  )
  
  curve_km <- build_curve_table(
    fit = km_transition,
    dataset_label = dataset_label,
    file_prefix = file_prefix,
    curve_type = "km_transition",
    max_followup_days = max_fu
  )
  
  curve_rkm <- build_curve_table(
    fit = reverse_km,
    dataset_label = dataset_label,
    file_prefix = file_prefix,
    curve_type = "reverse_km_followup",
    max_followup_days = max_fu
  )
  
  timepoints_all <- build_timepoint_table(
    fit = km_transition,
    times_days = eval_times_days,
    dataset_label = dataset_label,
    file_prefix = file_prefix,
    max_followup_days = max_fu
  )
  
  timepoints_primary <- timepoints_all |>
    dplyr::filter(within_observed_followup)
  
  summary_row <- build_summary_row(
    dat = dat,
    km_fit = km_transition,
    reverse_km = reverse_km,
    dataset_label = dataset_label,
    file_prefix = file_prefix
  )
  
  bundle <- list(
    meta = list(
      dataset_label = dataset_label,
      file_prefix = file_prefix,
      created_at = Sys.time(),
      eval_years = eval_years,
      eval_times_days = eval_times_days,
      key_definition = "site_id = paste(site, id, sep = '_')",
      status_definition = c(
        "0 = right_censoring",
        "1 = transition",
        "2 = remission"
      ),
      status_factor_levels = STATUS_LEVELS,
      primary_endpoint_definition = list(
        event = "transition (status_num == 1)",
        censoring = "right_censoring (0) + remission (2)"
      ),
      one_minus_km_definition = "1-KM under remission-as-censoring assumption; descriptive only; not CIF",
      reverse_km_definition = "Reverse-KM of potential follow-up time with event = status_num != 1; descriptive only; not a formal sufficient-follow-up test"
    ),
    step1_input = dat |>
      dplyr::select(dplyr::any_of(c(
        "site_id", "id", "site",
        "sex_num", "sex_fact",
        "date_birth", "date_entry",
        "age_int", "age_exact_entry", "age_exact_followup",
        "days_followup", "status_num", "status"
      ))),
    tables = list(
      summary = summary_row,
      curves = dplyr::bind_rows(curve_km, curve_rkm),
      timepoints_primary = timepoints_primary,
      timepoints_all = timepoints_all
    ),
    fit = list(
      km_transition = km_transition,
      reverse_km = reverse_km
    )
  )
  
  list(
    summary = summary_row,
    curves = dplyr::bind_rows(curve_km, curve_rkm),
    timepoints_primary = timepoints_primary,
    timepoints_all = timepoints_all,
    bundle = bundle
  )
}

## 🟠 Implement: plot rendering from curve data ===============================
render_curve_panel <- function(df, y_col, lower_col, upper_col, y_label, title_txt, subtitle_txt, out_png, show_censor_marks) {
  max_followup_days_vals <- unique(as.numeric(df$max_followup_days))
  max_followup_days_vals <- max_followup_days_vals[is.finite(max_followup_days_vals)]
  
  if (length(max_followup_days_vals) != 1L) {
    stop("max_followup_days must be unique within a plotting data frame.")
  }
  
  max_followup_days_scalar <- max_followup_days_vals[1]
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = time_days, y = .data[[y_col]])) +
    ggplot2::geom_step(linewidth = 0.9)
  
  if (!all(is.na(df[[lower_col]])) && !all(is.na(df[[upper_col]]))) {
    p <- p +
      ggplot2::geom_step(ggplot2::aes(y = .data[[lower_col]]), linetype = 2, alpha = 0.7) +
      ggplot2::geom_step(ggplot2::aes(y = .data[[upper_col]]), linetype = 2, alpha = 0.7)
  }
  
  if (isTRUE(show_censor_marks)) {
    censor_df <- df |>
      dplyr::filter(n_censor > 0)
    
    if (nrow(censor_df) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = censor_df,
          ggplot2::aes(x = time_days, y = .data[[y_col]]),
          inherit.aes = FALSE,
          shape = 3,
          size = 1.5
        )
    }
  }
  
  p <- p +
    ggplot2::geom_vline(
      xintercept = max_followup_days_scalar,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    ggplot2::annotate(
      "text",
      x = max_followup_days_scalar,
      y = 0.05,
      label = paste0("End of observed follow-up\n", round(max_followup_days_scalar, 1), " days"),
      hjust = 1.05,
      vjust = 0,
      size = 3
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Days since entry",
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 12)
  
  ggplot2::ggsave(
    filename = out_png,
    plot = p,
    width = plot_width_in,
    height = plot_height_in,
    dpi = plot_dpi
  )
}

render_dataset_plots <- function(curve_source_df, dataset_label_value, file_prefix_value, export_dir) {
  dataset_label_scalar <- as_scalar_character(dataset_label_value, "dataset_label_value")
  file_prefix_scalar <- as_scalar_character(file_prefix_value, "file_prefix_value")
  
  df_km <- curve_source_df |>
    dplyr::filter(
      dataset_label == dataset_label_scalar,
      file_prefix == file_prefix_scalar,
      curve_type == "km_transition"
    )
  
  df_rkm <- curve_source_df |>
    dplyr::filter(
      dataset_label == dataset_label_scalar,
      file_prefix == file_prefix_scalar,
      curve_type == "reverse_km_followup"
    )
  
  if (nrow(df_km) == 0L || nrow(df_rkm) == 0L) {
    stop("Plotting data frame is empty for dataset: ", dataset_label_scalar)
  }
  
  render_curve_panel(
    df = df_km,
    y_col = "surv",
    lower_col = "lower_95",
    upper_col = "upper_95",
    y_label = "Transition-free survival probability",
    title_txt = paste0("[", dataset_label_scalar, "] Transition-free survival (KM)"),
    subtitle_txt = "Primary endpoint: transition; remission treated as censoring",
    out_png = file.path(export_dir, paste0(file_prefix_scalar, "_step1_survival_transitionfree.png")),
    show_censor_marks = TRUE
  )
  
  render_curve_panel(
    df = df_km,
    y_col = "one_minus_km_transition",
    lower_col = "one_minus_km_transition_lcl95",
    upper_col = "one_minus_km_transition_ucl95",
    y_label = "1 - KM for transition\n(remission treated as censoring; not CIF)",
    title_txt = paste0("[", dataset_label_scalar, "] 1 - KM for transition"),
    subtitle_txt = "Descriptive only; this is not a CIF",
    out_png = file.path(export_dir, paste0(file_prefix_scalar, "_step1_one_minus_km_transition.png")),
    show_censor_marks = FALSE
  )
  
  render_curve_panel(
    df = df_rkm,
    y_col = "surv",
    lower_col = "lower_95",
    upper_col = "upper_95",
    y_label = "Probability potential follow-up exceeds t\n(reverse-KM; descriptive only)",
    title_txt = paste0("[", dataset_label_scalar, "] Reverse-KM follow-up distribution"),
    subtitle_txt = "Descriptive only; not a formal sufficient-follow-up test",
    out_png = file.path(export_dir, paste0(file_prefix_scalar, "_step1_reversekm_followup.png")),
    show_censor_marks = TRUE
  )
}

# 🔴 Execute: primary-endpoint Step 1 workflow ===============================

## 🟠 Read: merged cohort and standardize schema ===============================
data_all <- readr::read_csv(path_combined, show_col_types = FALSE)
data_all <- run_basic_qc(dat = data_all, strict_schema = strict_schema)

## 🟠 Compose: analysis datasets and planned filenames ===============================
dataset_specs <- make_dataset_specs(data_all)
target_files <- make_target_files(dataset_specs, export_dir)

if (length(target_files) > max_generated_files) {
  stop(
    "Planned output file count (", length(target_files),
    ") exceeds hard cap (", max_generated_files, ")."
  )
}

## 🟠 Remove: stale outputs for this Step 1 run ===============================
unlink(target_files[file.exists(target_files)], force = TRUE)

## 🟠 Fit: KM and reverse-KM bundles per dataset ===============================
result_list <- lapply(dataset_specs, function(spec) {
  fit_step1_dataset(
    dat = spec$data,
    dataset_label = spec$dataset_label,
    file_prefix = spec$file_prefix,
    eval_years = eval_years,
    eval_times_days = eval_times_days
  )
})
names(result_list) <- names(dataset_specs)

## 🟠 Bind: combined CSV source tables across datasets ===============================
summary_df <- dplyr::bind_rows(lapply(result_list, `[[`, "summary")) |>
  dplyr::arrange(dataset_label)

curve_source_df <- dplyr::bind_rows(lapply(result_list, `[[`, "curves")) |>
  dplyr::arrange(dataset_label, curve_type, time_days)

timepoints_primary_df <- dplyr::bind_rows(lapply(result_list, `[[`, "timepoints_primary")) |>
  dplyr::arrange(dataset_label, time_days)

timepoints_all_df <- dplyr::bind_rows(lapply(result_list, `[[`, "timepoints_all")) |>
  dplyr::arrange(dataset_label, time_days)

## 🟠 Write: combined CSV source-of-truth files ===============================
write_csv_safe(
  summary_df,
  file.path(export_dir, "step1_dataset_summary.csv")
)

write_csv_safe(
  curve_source_df,
  file.path(export_dir, "step1_curve_source.csv")
)

write_csv_safe(
  timepoints_primary_df,
  file.path(export_dir, "step1_km_timepoints_primary.csv")
)

write_csv_safe(
  timepoints_all_df,
  file.path(export_dir, "step1_km_timepoints_all.csv")
)

## 🟠 Save: per-dataset RDS fit bundles ===============================
for (spec_name in names(dataset_specs)) {
  spec <- dataset_specs[[spec_name]]
  bundle <- result_list[[spec_name]]$bundle
  
  saveRDS(
    bundle,
    file.path(export_dir, paste0(spec$file_prefix, "_step1_fit_bundle.rds"))
  )
}

## 🟠 Render: PNG panels from exported curve table ===============================
for (spec in dataset_specs) {
  render_dataset_plots(
    curve_source_df = curve_source_df,
    dataset_label_value = spec$dataset_label,
    file_prefix_value = spec$file_prefix,
    export_dir = export_dir
  )
}

## 🟠 Signal: completion status to console ===============================
message("Step 1 finished successfully.")
message("Generated files: ", length(target_files), " / ", max_generated_files)