# 🔴 Configure: Step2 paths, subgroup plan, and file budget ===============================
step1_dir <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦1.Step1_KM/attachments"
export_dir <- "/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.분석 방법 및 결과/🟦2.Step2_follow-up maturity, KM tail stability/attachments"

dataset_prefixes <- c("merged", "pnu", "snu")
subgroup_vars_merged <- c("sex_fact")
include_merged_subgroups <- TRUE

plot_width_in <- 11
plot_height_in <- 8
plot_dpi <- 300
max_generated_files <- 20L

upper_time_shift_frac <- 1e-8
default_eval_years <- 1:10
STATUS_LEVELS <- c("right_censoring", "remission", "transition")
EXPECTED_SITES <- c("PNU", "SNU")

# 🔴 Attach: packages and create export directory ===============================
required_pkgs <- c("survival", "readr", "dplyr", "tibble", "ggplot2")

to_install <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(to_install) > 0L) {
  install.packages(to_install)
}

invisible(lapply(required_pkgs, library, character.only = TRUE))
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# 🔴 Define: scalar helpers and schema guards ===============================
## 🟠 Implement: basic utility operators and value accessors ===============================
write_csv_safe <- function(x, path) {
  readr::write_csv(x, path, na = "")
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

safe_num_name <- function(x, name) {
  if (is.null(x)) return(NA_real_)
  if (is.data.frame(x)) {
    if (name %in% names(x)) return(as.numeric(x[[name]][1]))
    return(NA_real_)
  }
  if (is.matrix(x)) {
    if (!is.null(colnames(x)) && name %in% colnames(x)) return(as.numeric(x[1, name]))
    return(NA_real_)
  }
  if (!is.null(names(x)) && name %in% names(x)) {
    return(as.numeric(unname(x[[name]])))
  }
  NA_real_
}

safe_chr_name <- function(x, name) {
  if (is.null(x)) return(NA_character_)
  if (is.data.frame(x)) {
    if (name %in% names(x)) return(as.character(x[[name]][1]))
    return(NA_character_)
  }
  if (is.matrix(x)) {
    if (!is.null(colnames(x)) && name %in% colnames(x)) return(as.character(x[1, name]))
    return(NA_character_)
  }
  if (!is.null(names(x)) && name %in% names(x)) {
    return(as.character(unname(x[[name]])))
  }
  NA_character_
}

fill_na_with <- function(x, fill) {
  out <- x
  out[is.na(out)] <- fill[is.na(out)]
  out
}

sanitize_label <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x <- tolower(x)
  ifelse(nchar(x) == 0, "missing", x)
}

status_num_to_label <- function(x) {
  dplyr::case_when(
    x == 0L ~ "right_censoring",
    x == 1L ~ "transition",
    x == 2L ~ "remission",
    TRUE    ~ NA_character_
  )
}

trapz_base <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

survfit_surv_at_times <- function(fit, times) {
  s <- summary(fit, times = times, extend = TRUE)
  out <- s$surv
  if (length(out) == 0L) return(rep(1, length(times)))
  out
}

## 🟠 Implement: Step0-consistent input standardization ===============================
standardize_primary_input <- function(dat) {
  required_cols <- c("id", "site", "days_followup", "status_num")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns in Step1 input: ", paste(missing_cols, collapse = ", "))
  }
  
  dat <- dat |>
    dplyr::mutate(
      id = trimws(as.character(id)),
      site = toupper(trimws(as.character(site))),
      days_followup = as.numeric(days_followup),
      status_num = as.integer(status_num),
      site_id = if ("site_id" %in% names(dat)) trimws(as.character(site_id)) else paste(site, id, sep = "_")
    )
  
  if (any(is.na(dat$id) | dat$id == "")) stop("id contains missing or blank values.")
  if (any(is.na(dat$site) | dat$site == "")) stop("site contains missing or blank values.")
  if (any(is.na(dat$site_id) | dat$site_id == "")) stop("site_id contains missing or blank values.")
  if (any(is.na(dat$days_followup))) stop("days_followup contains NA.")
  if (any(dat$days_followup < 0)) stop("days_followup contains values < 0.")
  if (any(is.na(dat$status_num))) stop("status_num contains NA.")
  if (!all(dat$status_num %in% c(0L, 1L, 2L))) stop("status_num must be in {0,1,2}.")
  if (anyDuplicated(dat$site_id) > 0L) stop("Duplicated site_id detected in Step1 input.")
  
  dat <- dat |>
    dplyr::mutate(
      status = factor(status_num_to_label(status_num), levels = STATUS_LEVELS),
      time_primary = days_followup,
      event_primary = as.integer(status_num == 1L),
      censoring_primary = as.integer(status_num != 1L),
      status_primary = factor(
        dplyr::if_else(event_primary == 1L, "transition", "censored"),
        levels = c("censored", "transition")
      )
    )
  
  if ("sex_num" %in% names(dat)) {
    dat <- dat |>
      dplyr::mutate(
        sex_num = dplyr::if_else(sex_num %in% c(0, 1), as.integer(sex_num), NA_integer_)
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
  
  unexpected_sites <- setdiff(sort(unique(dat$site)), EXPECTED_SITES)
  if (length(unexpected_sites) > 0L) {
    warning(
      "Unexpected site values detected in Step2 input: ",
      paste(unexpected_sites, collapse = ", ")
    )
  }
  
  dat |>
    dplyr::arrange(time_primary, site_id)
}

## 🟠 Implement: source-of-truth curve exporters ===============================
survfit_curve_export_table <- function(fit, unit_meta, curve_type, max_followup_days) {
  s <- summary(fit, censored = TRUE)
  
  if (length(s$time) == 0L) {
    out <- tibble::tibble(
      unit_order = unit_meta$unit_order,
      unit_id = unit_meta$unit_id,
      unit_type = unit_meta$unit_type,
      dataset_label = unit_meta$dataset_label,
      parent_dataset = unit_meta$parent_dataset,
      group_var = unit_meta$group_var,
      group_level = unit_meta$group_level,
      display_label = unit_meta$display_label,
      file_prefix = unit_meta$file_prefix,
      fit_source = unit_meta$fit_source,
      curve_type = curve_type,
      time_days = 0,
      time_years = 0,
      n_risk = fit$n[1],
      n_event = 0,
      n_censor = 0,
      surv = 1,
      std_err = 0,
      lower_95 = 1,
      upper_95 = 1,
      max_followup_days = max_followup_days,
      max_followup_years = max_followup_days / 365.25
    )
  } else {
    surv_est <- s$surv
    se <- fill_na_with(s$std.err, rep(0, length(surv_est)))
    lower <- fill_na_with(s$lower, surv_est)
    upper <- fill_na_with(s$upper, surv_est)
    
    out <- tibble::tibble(
      unit_order = unit_meta$unit_order,
      unit_id = unit_meta$unit_id,
      unit_type = unit_meta$unit_type,
      dataset_label = unit_meta$dataset_label,
      parent_dataset = unit_meta$parent_dataset,
      group_var = unit_meta$group_var,
      group_level = unit_meta$group_level,
      display_label = unit_meta$display_label,
      file_prefix = unit_meta$file_prefix,
      fit_source = unit_meta$fit_source,
      curve_type = curve_type,
      time_days = c(0, s$time),
      time_years = c(0, s$time / 365.25),
      n_risk = c(fit$n[1], s$n.risk),
      n_event = c(0, s$n.event),
      n_censor = c(0, s$n.censor),
      surv = c(1, surv_est),
      std_err = c(0, se),
      lower_95 = c(1, lower),
      upper_95 = c(1, upper),
      max_followup_days = max_followup_days,
      max_followup_years = max_followup_days / 365.25
    )
  }
  
  out
}

# 🔴 Load: Step1 summary file and fit bundles ===============================
## 🟠 Read: combined Step1 CSV and per-dataset RDS objects ===============================
load_step1_assets <- function(step1_dir, dataset_prefixes) {
  summary_path <- file.path(step1_dir, "step1_dataset_summary.csv")
  bundle_paths <- file.path(step1_dir, paste0(dataset_prefixes, "_step1_fit_bundle.rds"))
  names(bundle_paths) <- dataset_prefixes
  
  missing_paths <- c(
    if (!file.exists(summary_path)) "step1_dataset_summary.csv" else character(0),
    names(bundle_paths)[!file.exists(bundle_paths)]
  )
  
  if (length(missing_paths) > 0L) {
    stop("Missing Step1 outputs: ", paste(missing_paths, collapse = ", "))
  }
  
  bundles <- lapply(bundle_paths, readRDS)
  names(bundles) <- dataset_prefixes
  
  list(
    summary_path = summary_path,
    bundle_paths = bundle_paths,
    summary_df = readr::read_csv(summary_path, show_col_types = FALSE),
    bundles = bundles
  )
}

## 🟠 Check: Step1 bundle consistency against Step1 summary CSV ===============================
check_step1_bundle <- function(bundle, step1_summary_df) {
  has_meta <- !is.null(bundle$meta)
  dataset_label <- if (has_meta) bundle$meta$dataset_label %||% NA_character_ else NA_character_
  file_prefix <- if (has_meta) bundle$meta$file_prefix %||% NA_character_ else NA_character_
  
  has_step1_input <- "step1_input" %in% names(bundle)
  has_fit <- "fit" %in% names(bundle)
  has_km_transition <- has_fit && "km_transition" %in% names(bundle$fit) && inherits(bundle$fit$km_transition, "survfit")
  has_reverse_km <- has_fit && "reverse_km" %in% names(bundle$fit) && inherits(bundle$fit$reverse_km, "survfit")
  has_summary_tbl <- "tables" %in% names(bundle) && "summary" %in% names(bundle$tables) && is.data.frame(bundle$tables$summary)
  
  input_n <- if (has_step1_input) nrow(bundle$step1_input) else NA_real_
  input_max_followup <- if (has_step1_input) max(as.numeric(bundle$step1_input$days_followup), na.rm = TRUE) else NA_real_
  
  bundle_summary_tbl <- if (has_summary_tbl) bundle$tables$summary else tibble::tibble()
  bundle_n <- safe_num_name(bundle_summary_tbl, "n")
  bundle_max_followup <- safe_num_name(bundle_summary_tbl, "max_followup_days")
  
  summary_match <- step1_summary_df |>
    dplyr::filter(file_prefix == !!file_prefix, dataset_label == !!dataset_label)
  
  summary_rows_found <- nrow(summary_match)
  summary_n <- safe_num_name(summary_match, "n")
  summary_max_followup <- safe_num_name(summary_match, "max_followup_days")
  
  endpoint_event <- if (has_meta && is.list(bundle$meta$primary_endpoint_definition)) {
    bundle$meta$primary_endpoint_definition$event %||% NA_character_
  } else {
    NA_character_
  }
  
  endpoint_censoring <- if (has_meta && is.list(bundle$meta$primary_endpoint_definition)) {
    bundle$meta$primary_endpoint_definition$censoring %||% NA_character_
  } else {
    NA_character_
  }
  
  step0_consistent <- identical(endpoint_event, "transition (status_num == 1)") &&
    identical(endpoint_censoring, "right_censoring (0) + remission (2)")
  
  qc_pass <- has_meta &&
    has_step1_input &&
    has_km_transition &&
    has_reverse_km &&
    has_summary_tbl &&
    summary_rows_found == 1L &&
    !is.na(input_n) &&
    !is.na(bundle_n) &&
    !is.na(summary_n) &&
    input_n == bundle_n &&
    input_n == summary_n &&
    !is.na(input_max_followup) &&
    !is.na(bundle_max_followup) &&
    !is.na(summary_max_followup) &&
    isTRUE(all.equal(input_max_followup, bundle_max_followup)) &&
    isTRUE(all.equal(input_max_followup, summary_max_followup)) &&
    step0_consistent
  
  tibble::tibble(
    qc_scope = "step1_bundle",
    unit_id = paste0("dataset__", file_prefix),
    dataset_label = dataset_label,
    parent_dataset = dataset_label,
    group_var = NA_character_,
    group_level = NA_character_,
    display_label = dataset_label,
    file_prefix = file_prefix,
    has_meta = has_meta,
    has_step1_input = has_step1_input,
    has_km_transition = has_km_transition,
    has_reverse_km = has_reverse_km,
    has_summary_tbl = has_summary_tbl,
    summary_rows_found = summary_rows_found,
    input_n = input_n,
    bundle_n = bundle_n,
    summary_n = summary_n,
    input_max_followup = input_max_followup,
    bundle_max_followup = bundle_max_followup,
    summary_max_followup = summary_max_followup,
    step0_consistent = step0_consistent,
    qc_pass = qc_pass,
    note = "Bundle QC against Step1 summary CSV and Step0 transition-only endpoint specification."
  )
}

# 🔴 Compose: analysis units from Step1 bundles ===============================
## 🟠 Build: dataset-level units by reusing Step1 bundle inputs and fits ===============================
make_dataset_unit_from_bundle <- function(bundle, unit_order) {
  dat <- standardize_primary_input(bundle$step1_input)
  
  list(
    unit_order = unit_order,
    unit_id = paste0("dataset__", bundle$meta$file_prefix),
    unit_type = "dataset",
    dataset_label = bundle$meta$dataset_label,
    parent_dataset = bundle$meta$dataset_label,
    group_var = NA_character_,
    group_level = NA_character_,
    display_label = bundle$meta$dataset_label,
    file_prefix = bundle$meta$file_prefix,
    fit_source = "reused_step1_fit",
    data = dat,
    prefit = bundle$fit
  )
}

## 🟠 Build: optional merged-cohort subgroup units for key strata ===============================
make_optional_subgroup_units <- function(merged_bundle, subgroup_vars, start_order) {
  units <- list()
  if (!length(subgroup_vars)) return(units)
  
  dat <- standardize_primary_input(merged_bundle$step1_input)
  order_counter <- start_order
  
  for (var_name in subgroup_vars) {
    if (!(var_name %in% names(dat))) next
    
    values <- as.character(dat[[var_name]])
    values[trimws(values) == ""] <- NA_character_
    levels_nonmissing <- unique(values[!is.na(values)])
    
    if (length(levels_nonmissing) == 0L) next
    
    for (lvl in levels_nonmissing) {
      subdat <- dat[values == lvl, , drop = FALSE]
      if (nrow(subdat) == 0L) next
      
      file_suffix <- paste0(sanitize_label(var_name), "_", sanitize_label(lvl))
      
      unit <- list(
        unit_order = order_counter,
        unit_id = paste("subgroup", "merged", sanitize_label(var_name), sanitize_label(lvl), sep = "__"),
        unit_type = "subgroup",
        dataset_label = "MERGED",
        parent_dataset = "MERGED",
        group_var = var_name,
        group_level = lvl,
        display_label = paste0("MERGED | ", var_name, "=", lvl),
        file_prefix = paste0("merged_", file_suffix),
        fit_source = "step2_refit_subgroup",
        data = subdat,
        prefit = NULL
      )
      
      units[[unit$unit_id]] <- unit
      order_counter <- order_counter + 1L
    }
  }
  
  units
}

## 🟠 Combine: dataset-level and subgroup-level analysis units ===============================
make_analysis_units <- function(step1_assets, include_merged_subgroups, subgroup_vars_merged) {
  units <- list()
  order_counter <- 1L
  
  for (prefix in names(step1_assets$bundles)) {
    units[[paste0("dataset__", prefix)]] <- make_dataset_unit_from_bundle(
      bundle = step1_assets$bundles[[prefix]],
      unit_order = order_counter
    )
    order_counter <- order_counter + 1L
  }
  
  if (isTRUE(include_merged_subgroups) && "merged" %in% names(step1_assets$bundles)) {
    subgroup_units <- make_optional_subgroup_units(
      merged_bundle = step1_assets$bundles[["merged"]],
      subgroup_vars = subgroup_vars_merged,
      start_order = order_counter
    )
    if (length(subgroup_units) > 0L) {
      units <- c(units, subgroup_units)
    }
  }
  
  unit_ids <- vapply(units, function(x) x$unit_id, character(1))
  if (anyDuplicated(unit_ids) > 0L) {
    stop("Duplicated unit_id detected in Step2 analysis units.")
  }
  
  units
}

# 🔴 Derive: descriptive summaries and KM stability limits ===============================
## 🟠 Compute: primary-tail summaries from KM and raw endpoint coding ===============================
compute_tail_summary <- function(dat, km_fit) {
  last_followup_days <- max(dat$time_primary, na.rm = TRUE)
  
  if (!any(dat$event_primary == 1L, na.rm = TRUE)) {
    return(
      tibble::tibble(
        last_transition_days = NA_real_,
        last_followup_days = last_followup_days,
        plateau_length_days = NA_real_,
        n_risk_at_last_transition = NA_real_,
        km_at_last_transition = NA_real_,
        one_more_event_drop_pct = NA_real_,
        n_right_censor_after_last_transition = NA_real_,
        n_remission_after_last_transition = NA_real_,
        n_primary_censor_after_last_transition = NA_real_
      )
    )
  }
  
  last_transition_days <- max(dat$time_primary[dat$event_primary == 1L], na.rm = TRUE)
  km_s <- summary(km_fit)
  
  idx <- max(which(km_s$time <= last_transition_days))
  n_risk_last <- km_s$n.risk[idx]
  surv_last <- km_s$surv[idx]
  
  tibble::tibble(
    last_transition_days = last_transition_days,
    last_followup_days = last_followup_days,
    plateau_length_days = last_followup_days - last_transition_days,
    n_risk_at_last_transition = n_risk_last,
    km_at_last_transition = surv_last,
    one_more_event_drop_pct = ifelse(
      is.na(n_risk_last) || n_risk_last <= 0,
      NA_real_,
      100 * surv_last / n_risk_last
    ),
    n_right_censor_after_last_transition = sum(dat$status_num == 0L & dat$time_primary > last_transition_days, na.rm = TRUE),
    n_remission_after_last_transition = sum(dat$status_num == 2L & dat$time_primary > last_transition_days, na.rm = TRUE),
    n_primary_censor_after_last_transition = sum(dat$status_num != 1L & dat$time_primary > last_transition_days, na.rm = TRUE)
  )
}

## 🟠 Compute: reverse-KM summaries for descriptive follow-up only ===============================
compute_reversekm_summary <- function(reverse_km) {
  tbl <- summary(reverse_km)$table
  
  tibble::tibble(
    reversekm_median_days = safe_num_name(tbl, "median"),
    reversekm_median_lcl95_days = safe_num_name(tbl, "0.95LCL"),
    reversekm_median_ucl95_days = safe_num_name(tbl, "0.95UCL"),
    reversekm_rmean_days = safe_num_name(tbl, "rmean"),
    reversekm_se_rmean_days = safe_num_name(tbl, "se(rmean)")
  )
}

## 🟠 Prepare: Betensky-style upper and lower limit datasets ===============================
make_betensky_limit_data <- function(dat, upper_time_shift_frac = 1e-8) {
  dat <- dat |>
    dplyr::transmute(
      site_id = as.character(site_id),
      time_primary = as.numeric(time_primary),
      event_primary = as.integer(event_primary)
    ) |>
    dplyr::arrange(time_primary)
  
  if (!any(dat$event_primary == 1L, na.rm = TRUE)) {
    return(list(
      last_transition_days = NA_real_,
      upper = NULL,
      lower = NULL,
      note = "No transition event: Betensky-style limits not available."
    ))
  }
  
  last_transition_days <- max(dat$time_primary[dat$event_primary == 1L], na.rm = TRUE)
  tiny_shift <- max(last_transition_days * upper_time_shift_frac, .Machine$double.eps)
  event_times <- sort(unique(dat$time_primary[dat$event_primary == 1L]))
  
  upper_dat <- dat
  idx_cens <- which(upper_dat$event_primary == 0L)
  if (length(idx_cens) > 0L) {
    upper_dat$time_primary[idx_cens] <- pmax(upper_dat$time_primary[idx_cens], last_transition_days + tiny_shift)
    upper_dat$event_primary[idx_cens] <- 0L
  }
  
  lower_dat <- dat
  if (length(idx_cens) > 0L) {
    next_event_time <- function(x, event_times) {
      cand <- event_times[event_times > x]
      if (length(cand) == 0L) NA_real_ else cand[1]
    }
    
    next_times <- vapply(
      lower_dat$time_primary[idx_cens],
      next_event_time,
      numeric(1),
      event_times = event_times
    )
    
    lower_dat$time_primary[idx_cens] <- ifelse(
      is.na(next_times),
      lower_dat$time_primary[idx_cens],
      next_times
    )
    lower_dat$event_primary[idx_cens] <- 1L
  }
  
  list(
    last_transition_days = last_transition_days,
    upper = upper_dat,
    lower = lower_dat,
    note = paste(
      "Betensky-style upper limit with exact post-last-event censoring shift;",
      "lower limit uses next-event recoding and falls back to censoring-time event recoding at the end of the tail (approximation)."
    )
  )
}

## 🟠 Fit: Betensky-style KM stability curves and summary metrics ===============================
fit_betensky_limits <- function(dat, km_fit, unit_meta, upper_time_shift_frac = 1e-8) {
  limdat <- make_betensky_limit_data(dat, upper_time_shift_frac = upper_time_shift_frac)
  max_followup_days <- max(dat$time_primary, na.rm = TRUE)
  
  if (is.null(limdat$upper) || is.null(limdat$lower)) {
    metrics_tbl <- tibble::tibble(
      unit_order = unit_meta$unit_order,
      unit_id = unit_meta$unit_id,
      unit_type = unit_meta$unit_type,
      dataset_label = unit_meta$dataset_label,
      parent_dataset = unit_meta$parent_dataset,
      group_var = unit_meta$group_var,
      group_level = unit_meta$group_level,
      display_label = unit_meta$display_label,
      file_prefix = unit_meta$file_prefix,
      fit_source = unit_meta$fit_source,
      last_transition_days = NA_real_,
      betensky_normalized_auc_to_last_transition = NA_real_,
      betensky_max_upper_minus_lower_to_last_transition = NA_real_,
      betensky_max_upper_minus_km_to_last_transition = NA_real_,
      betensky_max_km_minus_lower_to_last_transition = NA_real_,
      betensky_note = limdat$note
    )
    
    curve_tbl <- tibble::tibble(
      unit_order = unit_meta$unit_order,
      unit_id = unit_meta$unit_id,
      unit_type = unit_meta$unit_type,
      dataset_label = unit_meta$dataset_label,
      parent_dataset = unit_meta$parent_dataset,
      group_var = unit_meta$group_var,
      group_level = unit_meta$group_level,
      display_label = unit_meta$display_label,
      file_prefix = unit_meta$file_prefix,
      fit_source = unit_meta$fit_source,
      last_transition_days = NA_real_,
      time_days = NA_real_,
      time_years = NA_real_,
      km = NA_real_,
      upper = NA_real_,
      lower = NA_real_,
      upper_minus_lower = NA_real_,
      upper_minus_km = NA_real_,
      km_minus_lower = NA_real_,
      max_followup_days = max_followup_days,
      max_followup_years = max_followup_days / 365.25,
      betensky_note = limdat$note
    )
    
    return(list(metrics = metrics_tbl, curve_tbl = curve_tbl, fits = NULL))
  }
  
  upper_fit <- survival::survfit(
    survival::Surv(time_primary, event_primary) ~ 1,
    data = limdat$upper,
    conf.type = "none"
  )
  
  lower_fit <- survival::survfit(
    survival::Surv(time_primary, event_primary) ~ 1,
    data = limdat$lower,
    conf.type = "none"
  )
  
  all_times <- sort(unique(c(
    0,
    km_fit$time %||% numeric(0),
    upper_fit$time %||% numeric(0),
    lower_fit$time %||% numeric(0),
    max_followup_days
  )))
  
  surv_km <- survfit_surv_at_times(km_fit, all_times)
  surv_upper <- survfit_surv_at_times(upper_fit, all_times)
  surv_lower <- survfit_surv_at_times(lower_fit, all_times)
  
  curve_tbl <- tibble::tibble(
    unit_order = unit_meta$unit_order,
    unit_id = unit_meta$unit_id,
    unit_type = unit_meta$unit_type,
    dataset_label = unit_meta$dataset_label,
    parent_dataset = unit_meta$parent_dataset,
    group_var = unit_meta$group_var,
    group_level = unit_meta$group_level,
    display_label = unit_meta$display_label,
    file_prefix = unit_meta$file_prefix,
    fit_source = unit_meta$fit_source,
    last_transition_days = limdat$last_transition_days,
    time_days = all_times,
    time_years = all_times / 365.25,
    km = surv_km,
    upper = surv_upper,
    lower = surv_lower,
    upper_minus_lower = surv_upper - surv_lower,
    upper_minus_km = surv_upper - surv_km,
    km_minus_lower = surv_km - surv_lower,
    max_followup_days = max_followup_days,
    max_followup_years = max_followup_days / 365.25,
    betensky_note = limdat$note
  )
  
  metric_domain <- curve_tbl |>
    dplyr::filter(time_days <= limdat$last_transition_days)
  
  metrics_tbl <- tibble::tibble(
    unit_order = unit_meta$unit_order,
    unit_id = unit_meta$unit_id,
    unit_type = unit_meta$unit_type,
    dataset_label = unit_meta$dataset_label,
    parent_dataset = unit_meta$parent_dataset,
    group_var = unit_meta$group_var,
    group_level = unit_meta$group_level,
    display_label = unit_meta$display_label,
    file_prefix = unit_meta$file_prefix,
    fit_source = unit_meta$fit_source,
    last_transition_days = limdat$last_transition_days,
    betensky_normalized_auc_to_last_transition = if (limdat$last_transition_days > 0) {
      trapz_base(metric_domain$time_days, metric_domain$upper_minus_lower) / limdat$last_transition_days
    } else {
      NA_real_
    },
    betensky_max_upper_minus_lower_to_last_transition = if (nrow(metric_domain) > 0L) max(metric_domain$upper_minus_lower, na.rm = TRUE) else NA_real_,
    betensky_max_upper_minus_km_to_last_transition = if (nrow(metric_domain) > 0L) max(metric_domain$upper_minus_km, na.rm = TRUE) else NA_real_,
    betensky_max_km_minus_lower_to_last_transition = if (nrow(metric_domain) > 0L) max(metric_domain$km_minus_lower, na.rm = TRUE) else NA_real_,
    betensky_note = limdat$note
  )
  
  list(
    metrics = metrics_tbl,
    curve_tbl = curve_tbl,
    fits = list(km = km_fit, upper = upper_fit, lower = lower_fit)
  )
}

## 🟠 Assemble: unit-level descriptive follow-up maturity row ===============================
assemble_descriptive_row <- function(unit_meta, dat, reverse_km, km_fit, betensky_metrics) {
  reverse_tbl <- compute_reversekm_summary(reverse_km)
  tail_tbl <- compute_tail_summary(dat, km_fit)
  
  counts_tbl <- tibble::tibble(
    unit_order = unit_meta$unit_order,
    unit_id = unit_meta$unit_id,
    unit_type = unit_meta$unit_type,
    dataset_label = unit_meta$dataset_label,
    parent_dataset = unit_meta$parent_dataset,
    group_var = unit_meta$group_var,
    group_level = unit_meta$group_level,
    display_label = unit_meta$display_label,
    file_prefix = unit_meta$file_prefix,
    fit_source = unit_meta$fit_source,
    n = nrow(dat),
    n_right_censoring = sum(dat$status_num == 0L, na.rm = TRUE),
    n_transition = sum(dat$status_num == 1L, na.rm = TRUE),
    n_remission = sum(dat$status_num == 2L, na.rm = TRUE),
    n_primary_events = sum(dat$event_primary == 1L, na.rm = TRUE),
    n_primary_censored = sum(dat$event_primary == 0L, na.rm = TRUE),
    min_followup_days = min(dat$time_primary, na.rm = TRUE),
    median_observed_time_days_raw = stats::median(dat$time_primary, na.rm = TRUE),
    max_followup_days = max(dat$time_primary, na.rm = TRUE),
    max_followup_years = max(dat$time_primary, na.rm = TRUE) / 365.25,
    remission_treated_as_censoring = TRUE,
    primary_endpoint_note = "Transition-only right-censored endpoint: remission and right-censoring are both treated as censoring.",
    reversekm_note = "Reverse-KM is descriptive only and not a formal sufficient-follow-up test.",
    betensky_note_interpretation = "Betensky-style stability limits are descriptive only; the lower limit uses an approximation at the end of the tail."
  )
  
  counts_tbl |>
    dplyr::bind_cols(reverse_tbl, tail_tbl, betensky_metrics |>
                       dplyr::select(
                         betensky_normalized_auc_to_last_transition,
                         betensky_max_upper_minus_lower_to_last_transition,
                         betensky_max_upper_minus_km_to_last_transition,
                         betensky_max_km_minus_lower_to_last_transition,
                         betensky_note
                       ))
}

# 🔴 Execute: per-unit descriptive Step2 analyses ===============================
## 🟠 Fit: reuse Step1 survival objects or refit subgroup units ===============================
fit_unit_surv_objects <- function(unit) {
  dat <- unit$data
  
  reuse_ok <- !is.null(unit$prefit) &&
    "km_transition" %in% names(unit$prefit) &&
    "reverse_km" %in% names(unit$prefit) &&
    inherits(unit$prefit$km_transition, "survfit") &&
    inherits(unit$prefit$reverse_km, "survfit")
  
  if (isTRUE(reuse_ok)) {
    unit$fit_source <- "reused_step1_fit"
    return(list(
      km_transition = unit$prefit$km_transition,
      reverse_km = unit$prefit$reverse_km,
      fit_source = unit$fit_source
    ))
  }
  
  km_transition <- survival::survfit(
    survival::Surv(time_primary, event_primary) ~ 1,
    data = dat,
    conf.type = "log-log"
  )
  
  reverse_km <- survival::survfit(
    survival::Surv(time_primary, censoring_primary) ~ 1,
    data = dat,
    conf.type = "log-log"
  )
  
  unit$fit_source <- "step2_refit_subgroup"
  
  list(
    km_transition = km_transition,
    reverse_km = reverse_km,
    fit_source = unit$fit_source
  )
}

## 🟠 Check: unit-level descriptive input integrity before export ===============================
check_analysis_unit <- function(unit) {
  dat <- unit$data
  
  qc_pass <- nrow(dat) > 0L &&
    all(!is.na(dat$site_id)) &&
    anyDuplicated(dat$site_id) == 0L &&
    all(!is.na(dat$time_primary)) &&
    all(dat$time_primary >= 0) &&
    all(!is.na(dat$status_num)) &&
    all(dat$status_num %in% c(0L, 1L, 2L))
  
  tibble::tibble(
    qc_scope = "analysis_unit",
    unit_id = unit$unit_id,
    dataset_label = unit$dataset_label,
    parent_dataset = unit$parent_dataset,
    group_var = unit$group_var,
    group_level = unit$group_level,
    display_label = unit$display_label,
    file_prefix = unit$file_prefix,
    unit_order = unit$unit_order,
    n = nrow(dat),
    n_transition = sum(dat$event_primary == 1L, na.rm = TRUE),
    n_primary_censored = sum(dat$event_primary == 0L, na.rm = TRUE),
    max_followup_days = max(dat$time_primary, na.rm = TRUE),
    unique_site_ids = length(unique(dat$site_id)),
    qc_pass = qc_pass,
    note = "Analysis-unit QC under the Step0 transition-only endpoint definition."
  )
}

## 🟠 Run: single-unit Step2 descriptive workflow and collect reusable objects ===============================
run_step2_one_unit <- function(unit, upper_time_shift_frac = 1e-8) {
  qc_tbl <- check_analysis_unit(unit)
  if (!isTRUE(qc_tbl$qc_pass[1])) {
    stop("Analysis-unit QC failed for: ", unit$display_label)
  }
  
  fits <- fit_unit_surv_objects(unit)
  unit$fit_source <- fits$fit_source
  max_followup_days <- max(unit$data$time_primary, na.rm = TRUE)
  
  reverse_curve_tbl <- survfit_curve_export_table(
    fit = fits$reverse_km,
    unit_meta = unit,
    curve_type = "reverse_km_followup",
    max_followup_days = max_followup_days
  )
  
  betensky_obj <- fit_betensky_limits(
    dat = unit$data,
    km_fit = fits$km_transition,
    unit_meta = unit,
    upper_time_shift_frac = upper_time_shift_frac
  )
  
  descriptive_tbl <- assemble_descriptive_row(
    unit_meta = unit,
    dat = unit$data,
    reverse_km = fits$reverse_km,
    km_fit = fits$km_transition,
    betensky_metrics = betensky_obj$metrics
  )
  
  list(
    qc = qc_tbl,
    descriptive = descriptive_tbl,
    reverse_curve = reverse_curve_tbl,
    betensky_curve = betensky_obj$curve_tbl,
    fits = list(
      km_transition = fits$km_transition,
      reverse_km = fits$reverse_km,
      betensky = betensky_obj$fits
    ),
    data = unit$data
  )
}

# 🔴 Render: faceted Step2 plots from exported CSV tables only ===============================
## 🟠 Draw: reverse-KM follow-up panels from combined CSV ===============================
render_reversekm_plot_from_csv <- function(curve_csv, out_png) {
  curve_df <- readr::read_csv(curve_csv, show_col_types = FALSE)
  
  if (nrow(curve_df) == 0L) {
    ggplot2::ggsave(
      filename = out_png,
      plot = ggplot2::ggplot() + ggplot2::theme_void(),
      width = plot_width_in,
      height = plot_height_in,
      dpi = plot_dpi
    )
    return(invisible(NULL))
  }
  
  label_order <- curve_df |>
    dplyr::distinct(unit_order, display_label) |>
    dplyr::arrange(unit_order) |>
    dplyr::pull(display_label)
  
  curve_df <- curve_df |>
    dplyr::mutate(
      display_label = factor(display_label, levels = label_order)
    )
  
  vline_df <- curve_df |>
    dplyr::distinct(unit_id, unit_order, display_label, max_followup_days) |>
    dplyr::arrange(unit_order)
  
  p <- ggplot2::ggplot(curve_df, ggplot2::aes(x = time_days, y = surv)) +
    ggplot2::geom_step(linewidth = 0.85) +
    ggplot2::geom_step(ggplot2::aes(y = lower_95), linetype = 2, alpha = 0.7) +
    ggplot2::geom_step(ggplot2::aes(y = upper_95), linetype = 2, alpha = 0.7) +
    ggplot2::geom_vline(
      data = vline_df,
      ggplot2::aes(xintercept = max_followup_days),
      inherit.aes = FALSE,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::facet_wrap(~ display_label, scales = "free_x", ncol = 2) +
    ggplot2::labs(
      title = "Step 2 reverse-KM follow-up distribution",
      subtitle = "Descriptive only; dashed vertical line marks the maximum observed follow-up in each panel",
      x = "Days since entry",
      y = "Probability potential follow-up exceeds t"
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

## 🟠 Draw: Betensky-style KM/upper/lower overlays from combined CSV ===============================
render_betensky_overlay_plot_from_csv <- function(curve_csv, out_png) {
  curve_df <- readr::read_csv(curve_csv, show_col_types = FALSE) |>
    dplyr::filter(!is.na(km))
  
  if (nrow(curve_df) == 0L) {
    ggplot2::ggsave(
      filename = out_png,
      plot = ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No transition events: no Betensky-style curves available.") +
        ggplot2::theme_void(),
      width = plot_width_in,
      height = plot_height_in,
      dpi = plot_dpi
    )
    return(invisible(NULL))
  }
  
  label_order <- curve_df |>
    dplyr::distinct(unit_order, display_label) |>
    dplyr::arrange(unit_order) |>
    dplyr::pull(display_label)
  
  curve_df <- curve_df |>
    dplyr::mutate(
      display_label = factor(display_label, levels = label_order)
    )
  
  long_df <- dplyr::bind_rows(
    curve_df |>
      dplyr::transmute(unit_order, display_label, time_days, curve = "KM", surv = km),
    curve_df |>
      dplyr::transmute(unit_order, display_label, time_days, curve = "Upper limit", surv = upper),
    curve_df |>
      dplyr::transmute(unit_order, display_label, time_days, curve = "Lower limit", surv = lower)
  )
  
  last_event_df <- curve_df |>
    dplyr::filter(!is.na(last_transition_days)) |>
    dplyr::distinct(unit_order, display_label, last_transition_days)
  
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = time_days, y = surv, linetype = curve)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::geom_vline(
      data = last_event_df,
      ggplot2::aes(xintercept = last_transition_days),
      inherit.aes = FALSE,
      linetype = "dotted",
      linewidth = 0.5
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::facet_wrap(~ display_label, scales = "free_x", ncol = 2) +
    ggplot2::labs(
      title = "Step 2 KM tail stability with Betensky-style limits",
      subtitle = "Descriptive only; dotted vertical line marks the last observed transition time in each panel",
      x = "Days since entry",
      y = "Transition-free survival probability",
      linetype = NULL
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

## 🟠 Draw: Betensky-style difference curves up to last transition time ===============================
render_betensky_difference_plot_from_csv <- function(curve_csv, out_png) {
  curve_df <- readr::read_csv(curve_csv, show_col_types = FALSE) |>
    dplyr::filter(!is.na(km), !is.na(last_transition_days), time_days <= last_transition_days)
  
  if (nrow(curve_df) == 0L) {
    ggplot2::ggsave(
      filename = out_png,
      plot = ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No transition events: no difference curves available.") +
        ggplot2::theme_void(),
      width = plot_width_in,
      height = plot_height_in,
      dpi = plot_dpi
    )
    return(invisible(NULL))
  }
  
  label_order <- curve_df |>
    dplyr::distinct(unit_order, display_label) |>
    dplyr::arrange(unit_order) |>
    dplyr::pull(display_label)
  
  curve_df <- curve_df |>
    dplyr::mutate(
      display_label = factor(display_label, levels = label_order)
    )
  
  long_df <- dplyr::bind_rows(
    curve_df |>
      dplyr::transmute(unit_order, display_label, time_days, measure = "Upper - Lower", value = upper_minus_lower),
    curve_df |>
      dplyr::transmute(unit_order, display_label, time_days, measure = "Upper - KM", value = upper_minus_km),
    curve_df |>
      dplyr::transmute(unit_order, display_label, time_days, measure = "KM - Lower", value = km_minus_lower)
  )
  
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = time_days, y = value, linetype = measure)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::facet_wrap(~ display_label, scales = "free_x", ncol = 2) +
    ggplot2::labs(
      title = "Step 2 Betensky-style difference curves",
      subtitle = "Restricted to the time range up to the last observed transition in each panel",
      x = "Days since entry",
      y = "Probability difference",
      linetype = NULL
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

# 🔴 Plan: Step2 target files and stale-output cleanup ===============================
## 🟠 Declare: fixed output manifest under the file-budget rule ===============================
target_files <- file.path(export_dir, c(
  "step2_input_qc.csv",
  "step2_analysis_unit_registry.csv",
  "step2_followup_descriptive.csv",
  "step2_reversekm_curve.csv",
  "step2_km_stability_limits.csv",
  "step2_descriptive_bundle.rds",
  "step2_reversekm_followup.png",
  "step2_km_stability_overlay.png",
  "step2_km_stability_difference.png",
  "step2_manifest.csv"
))

if (length(target_files) > max_generated_files) {
  stop(
    "Planned Step2 output file count (", length(target_files),
    ") exceeds hard cap (", max_generated_files, ")."
  )
}

unlink(target_files[file.exists(target_files)], force = TRUE)

# 🔴 Run: Step2 descriptive-only workflow under Step0 endpoint lock ===============================
## 🟠 Load: Step1 inputs and verify transition-only bundle compatibility ===============================
step1_assets <- load_step1_assets(
  step1_dir = step1_dir,
  dataset_prefixes = dataset_prefixes
)

bundle_qc_df <- dplyr::bind_rows(lapply(step1_assets$bundles, check_step1_bundle, step1_summary_df = step1_assets$summary_df)) |>
  dplyr::arrange(dataset_label)

if (!all(bundle_qc_df$qc_pass)) {
  stop("At least one Step1 bundle failed QC. Inspect step2_input_qc.csv after fixing upstream outputs.")
}

## 🟠 Create: dataset-level and subgroup-level units for descriptive Step2 ===============================
analysis_units <- make_analysis_units(
  step1_assets = step1_assets,
  include_merged_subgroups = include_merged_subgroups,
  subgroup_vars_merged = subgroup_vars_merged
)

unit_results <- lapply(
  analysis_units,
  function(unit) {
    run_step2_one_unit(
      unit = unit,
      upper_time_shift_frac = upper_time_shift_frac
    )
  }
)
names(unit_results) <- names(analysis_units)

## 🟠 Bind: source-of-truth tables across all Step2 analysis units ===============================
unit_qc_df <- dplyr::bind_rows(lapply(unit_results, `[[`, "qc")) |>
  dplyr::arrange(unit_order)

input_qc_df <- dplyr::bind_rows(bundle_qc_df, unit_qc_df)

followup_descriptive_df <- dplyr::bind_rows(lapply(unit_results, `[[`, "descriptive")) |>
  dplyr::arrange(unit_order)

analysis_unit_registry_df <- followup_descriptive_df |>
  dplyr::select(
    unit_order,
    unit_id,
    unit_type,
    dataset_label,
    parent_dataset,
    group_var,
    group_level,
    display_label,
    file_prefix,
    fit_source,
    n,
    n_transition,
    n_primary_censored,
    max_followup_days
  ) |>
  dplyr::arrange(unit_order)

reversekm_curve_df <- dplyr::bind_rows(lapply(unit_results, `[[`, "reverse_curve")) |>
  dplyr::arrange(unit_order, time_days)

km_stability_limits_df <- dplyr::bind_rows(lapply(unit_results, `[[`, "betensky_curve")) |>
  dplyr::arrange(unit_order, time_days)

# 🔴 Export: Step2 CSV tables and reusable RDS bundle ===============================
## 🟠 Write: combined descriptive source-of-truth files ===============================
write_csv_safe(
  input_qc_df,
  file.path(export_dir, "step2_input_qc.csv")
)

write_csv_safe(
  analysis_unit_registry_df,
  file.path(export_dir, "step2_analysis_unit_registry.csv")
)

write_csv_safe(
  followup_descriptive_df,
  file.path(export_dir, "step2_followup_descriptive.csv")
)

write_csv_safe(
  reversekm_curve_df,
  file.path(export_dir, "step2_reversekm_curve.csv")
)

write_csv_safe(
  km_stability_limits_df,
  file.path(export_dir, "step2_km_stability_limits.csv")
)

## 🟠 Save: self-contained Step2 RDS for downstream survival summaries ===============================
step2_bundle <- list(
  meta = list(
    created_at = Sys.time(),
    step0_primary_endpoint_definition = list(
      event = "transition (status_num == 1)",
      censoring = "right_censoring (0) + remission (2)"
    ),
    step2_role = "Descriptive-only assessment of follow-up maturity and KM tail stability before any formal cure-model gate.",
    methods = list(
      reverse_km = "Descriptive reverse-KM summary of potential follow-up time",
      tail_summary = "Last transition, last follow-up, plateau length, n at risk, and censoring counts after the last transition",
      betensky_style_limits = "Betensky-style upper limit plus lower-limit approximation for KM tail stability"
    ),
    assumptions = list(
      remission_treated_as_censoring = TRUE,
      right_censoring_framework = TRUE,
      no_formal_gate_in_step2 = TRUE
    ),
    default_eval_years = default_eval_years,
    default_eval_times_days = default_eval_years * 365.25,
    source_step1_dir = step1_dir,
    source_step1_summary_csv = step1_assets$summary_path,
    source_step1_bundle_paths = step1_assets$bundle_paths
  ),
  step1_reference = list(
    summary_df = step1_assets$summary_df,
    bundle_qc = bundle_qc_df
  ),
  analysis_units = analysis_unit_registry_df,
  input_qc = input_qc_df,
  data = lapply(unit_results, `[[`, "data"),
  fit = lapply(unit_results, `[[`, "fits"),
  results = list(
    followup_descriptive = followup_descriptive_df,
    reversekm_curve = reversekm_curve_df,
    km_stability_limits = km_stability_limits_df
  )
)

saveRDS(
  step2_bundle,
  file.path(export_dir, "step2_descriptive_bundle.rds")
)

# 🔴 Render: PNG panels strictly from exported CSV tables ===============================
## 🟠 Generate: reverse-KM, stability overlay, and difference plots ===============================
render_reversekm_plot_from_csv(
  curve_csv = file.path(export_dir, "step2_reversekm_curve.csv"),
  out_png = file.path(export_dir, "step2_reversekm_followup.png")
)

render_betensky_overlay_plot_from_csv(
  curve_csv = file.path(export_dir, "step2_km_stability_limits.csv"),
  out_png = file.path(export_dir, "step2_km_stability_overlay.png")
)

render_betensky_difference_plot_from_csv(
  curve_csv = file.path(export_dir, "step2_km_stability_limits.csv"),
  out_png = file.path(export_dir, "step2_km_stability_difference.png")
)

# 🔴 Register: final manifest after all exports and plots ===============================
## 🟠 Write: file inventory with size and timestamp metadata ===============================
manifest_df <- tibble::tibble(
  file = basename(target_files[file.exists(target_files)])
) |>
  dplyr::mutate(
    full_path = file.path(export_dir, file),
    size_bytes = file.info(full_path)$size,
    modified_time = as.character(file.info(full_path)$mtime)
  ) |>
  dplyr::arrange(file)

write_csv_safe(
  manifest_df,
  file.path(export_dir, "step2_manifest.csv")
)

# 🔴 Signal: completion status to console ===============================
message("Step 2 finished successfully.")
message("Generated files: ", length(target_files), " / ", max_generated_files)