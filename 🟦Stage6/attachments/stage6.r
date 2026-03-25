# 🔴 Stage 6 cure appropriateness screening — spec-aligned SHARDED runtime (hotfix integrated) ===============================
# This script is designed so that the SAME FILE can be launched in multiple R windows concurrently.
# Each worker claims pending bootstrap shards from a shared registry, computes them, writes shard RDS files,
# and one worker performs the final reduce/export step after all shards are complete.

# 🔴 Configure: run-root, input paths, and shard controls ===============================
run_mode <- "full"  # one of: smoke, balanced, full

run_root_dir <- '/Users/ido/Library/CloudStorage/Dropbox/Data Analysis/STAGE6/stage6_run__2026-03-25_hsuShard10_xieShard50_auto_v1__000001'
export_path <- run_root_dir

# Optional: if TRUE, copy the final canonical Stage 6 outputs into the project attachments folder after reduction.
copy_final_exports_to_canonical <- FALSE
canonical_export_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage6/attachments'

# Use a fresh run_root_dir for each new run. Only set this to TRUE when you intentionally want to wipe the
# current run_root_dir state and restart the sharded run from scratch.
refresh_run_state <- FALSE
refresh_final_exports <- TRUE

# Worker-loop controls. The same script can be opened in 10 R windows at once.
worker_poll_seconds <- 5
max_worker_idle_cycles <- 720L
worker_nonblocking_reduce_lock_timeout_sec <- 1
worker_claim_lock_timeout_sec <- 600
worker_reduce_lock_timeout_sec <- 600

# Hotfix controls
task_stale_after_sec <- 3600L
failed_task_retry_after_sec <- 60L
max_task_attempts <- 3L
lock_stale_after_sec <- 21600L
task_heartbeat_every_reps <- 1L

# Shard sizes. With 10 simultaneous R windows, these settings usually work well for full mode.
xie_shard_size <- 50L
hsu_shard_size <- 10L

# Keep shard workers single-core when launching many R windows simultaneously.
parallel_cores <- 1L
allow_parallel_bootstrap <- FALSE
bootstrap_chunk_size <- 1L

# Input paths
data_path <- '/Volumes/ObsidianVault/Obsidian/☔️Papers_Writing(논문 쓰기)/📙Currently working/⬛조현병 베이지안 생존분석/🟧1.데이터 분석 New/🟦Stage1/attachments'
stage1_bundle_file <- file.path(data_path, "stage1_backbone_bundle.rds")
stage1_datasets_file <- file.path(data_path, "stage1_analysis_datasets.rds")
stage1_dataset_registry_file <- file.path(data_path, "stage1_dataset_registry.csv")
stage1_scaling_registry_file <- file.path(data_path, "stage1_scaling_registry.csv")
stage1_formula_registry_file <- file.path(data_path, "stage1_formula_registry.csv")
stage1_modeling_registry_file <- file.path(data_path, "stage1_modeling_registry.csv")
stage1_horizon_registry_file <- file.path(data_path, "stage1_horizon_registry.csv")
stage1_threshold_registry_file <- file.path(data_path, "stage1_threshold_registry.csv")
stage1_metadata_registry_file <- file.path(data_path, "stage1_metadata_registry.csv")

# Stage 6 screening constants
alpha_screening <- 0.05
bootstrap_seed <- 20260322L
receus_pi_threshold <- 0.025
receus_r_threshold <- 0.05
receus_pi_equivocal <- 0.010
receus_r_equivocal <- 0.10
hsu_min_events <- 5L
tiny_time <- 1e-8
save_fit_objects <- TRUE
max_log_param_abs <- 30

candidate_receus_families <- c("exponential", "weibull", "gamma", "loglogistic")
candidate_hsu_families <- c("weibull", "lognormal", "loglogistic")
hsu_familyset_method_name <- "hsu_supscore_familyset_approx"
hsu_family_set_name <- paste(candidate_hsu_families, collapse = "|")
stage6_specification_file <- "Stage6_Cure_Appropriateness_Screening_Spec_FINAL.md"
decision_rule_version <- "stage6_receus_aic_primary_hsu_familyset_xie_contradiction_v3"
screening_context <- "observational_cohort"

if (identical(run_mode, "smoke")) {
  bootstrap_reps_xie <- 19L
  bootstrap_reps_hsu <- 9L
  optim_maxit <- 400L
  optim_reltol <- 1e-8
  max_hsu_formula_count <- 2L
  xie_shard_size <- min(xie_shard_size, 10L)
  hsu_shard_size <- min(hsu_shard_size, 5L)
} else if (identical(run_mode, "balanced")) {
  bootstrap_reps_xie <- 99L
  bootstrap_reps_hsu <- 39L
  optim_maxit <- 1500L
  optim_reltol <- 1e-9
  max_hsu_formula_count <- 3L
  xie_shard_size <- min(xie_shard_size, 25L)
  hsu_shard_size <- min(hsu_shard_size, 5L)
} else if (identical(run_mode, "full")) {
  bootstrap_reps_xie <- 499L
  bootstrap_reps_hsu <- 199L
  optim_maxit <- 5000L
  optim_reltol <- 1e-10
  max_hsu_formula_count <- Inf
} else {
  stop("`run_mode` must be one of: smoke, balanced, full.", call. = FALSE)
}

# 🔴 Initialize: packages, directories, and shared run-state paths ===============================
required_packages <- c("survival", "dplyr", "readr", "tibble", "purrr", "stringr", "ggplot2", "scales")

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop(
    sprintf(
      "Install required packages before running Stage 6: %s",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 999)

dir.create(run_root_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

state_dir <- file.path(run_root_dir, "state")
observed_dir <- file.path(run_root_dir, "observed")
shards_dir <- file.path(run_root_dir, "shards")
xie_shards_dir <- file.path(shards_dir, "xie")
hsu_shards_dir <- file.path(shards_dir, "hsu")
merged_dir <- file.path(run_root_dir, "merged")
locks_dir <- file.path(run_root_dir, "_locks")

dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(observed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(xie_shards_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hsu_shards_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(locks_dir, recursive = TRUE, showWarnings = FALSE)

heartbeat_dir <- file.path(state_dir, "_heartbeats")
dir.create(heartbeat_dir, recursive = TRUE, showWarnings = FALSE)

progress_log_file <- file.path(run_root_dir, "stage6_progress_log.txt")
runtime_status_file <- file.path(run_root_dir, "stage6_runtime_status.csv")
shard_registry_file <- file.path(run_root_dir, "stage6_shard_registry.csv")
run_manifest_file <- file.path(run_root_dir, "stage6_run_manifest.csv")
state_bundle_file <- file.path(state_dir, "stage6_state_bundle.rds")
prepared_datasets_file <- file.path(state_dir, "stage6_prepared_datasets.rds")
state_initialized_file <- file.path(state_dir, "stage6_state_initialized.flag")
reduction_complete_file <- file.path(state_dir, "stage6_reduction_completed.flag")

init_lock_dir <- file.path(locks_dir, "init.lock")
registry_lock_dir <- file.path(locks_dir, "registry.lock")
reduce_lock_dir <- file.path(locks_dir, "reduce.lock")
runtime_lock_dir <- file.path(locks_dir, "runtime.lock")

worker_id <- paste0(
  Sys.info()[["nodename"]],
  "_pid",
  Sys.getpid(),
  "_",
  format(Sys.time(), "%H%M%S")
)

# 🔴 Define: generic runtime and I/O helpers ===============================
now_utc_chr <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

parse_timestamp_utc <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- NA_character_
  suppressWarnings(
    as.POSIXct(
      x,
      tz = "UTC",
      tryFormats = c(
        "%Y-%m-%dT%H:%M:%SZ",
        "%Y-%m-%d %H:%M:%S",
        "%Y/%m/%d %H:%M:%S"
      )
    )
  )
}

timestamp_age_sec <- function(x, default = Inf) {
  ts <- parse_timestamp_utc(x)
  if (length(ts) == 0L || is.na(ts[[1]])) {
    return(default)
  }
  as.numeric(difftime(Sys.time(), ts[[1]], units = "secs"))
}

as_integer_or_na <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

as_numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

coalesce_chr_scalar <- function(...) {
  xs <- list(...)
  for (ii in seq_along(xs)) {
    val <- xs[[ii]]
    if (length(val) == 0L || is.null(val)) {
      next
    }
    val <- as.character(val[[1]])
    if (!is.na(val) && nzchar(val)) {
      return(val)
    }
  }
  NA_character_
}

format_elapsed <- function(start_time) {
  elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (!is.finite(elapsed_sec)) {
    return("NA sec")
  }
  if (elapsed_sec < 60) {
    return(sprintf("%.1f sec", elapsed_sec))
  }
  if (elapsed_sec < 3600) {
    return(sprintf("%.1f min", elapsed_sec / 60))
  }
  sprintf("%.2f hr", elapsed_sec / 3600)
}

log_step <- function(..., level = "INFO") {
  line <- sprintf(
    "[%s] [%s] [%s] %s",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    level,
    worker_id,
    paste0(..., collapse = "")
  )
  message(line)
  cat(line, "\n", file = progress_log_file, append = TRUE)
  flush.console()
  invisible(line)
}

safe_read_csv_guess <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_read_csv_text <- function(path) {
  readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
}

lock_info_file <- function(lock_dir) {
  file.path(lock_dir, "lock_info.csv")
}

lock_timestamp_from_dir <- function(lock_dir) {
  if (file.exists(lock_info_file(lock_dir))) {
    info <- tryCatch(safe_read_csv_text(lock_info_file(lock_dir)), error = function(e) tibble::tibble())
    if (nrow(info) > 0L && "acquired_at_utc" %in% names(info)) {
      return(info$acquired_at_utc[[1]])
    }
  }
  info <- suppressWarnings(file.info(lock_dir))
  if (nrow(info) == 0L || is.na(info$mtime[[1]])) {
    return(NA_character_)
  }
  format(info$mtime[[1]], "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

write_lock_metadata <- function(lock_dir) {
  metadata <- tibble::tibble(
    worker_id = worker_id,
    pid = as.character(Sys.getpid()),
    acquired_at_utc = now_utc_chr()
  )
  temp_path <- tempfile(pattern = "lock_meta_", tmpdir = lock_dir, fileext = ".csv")
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  readr::write_csv(metadata, temp_path)
  file.rename(temp_path, lock_info_file(lock_dir))
  invisible(TRUE)
}

break_stale_lock_if_needed <- function(lock_dir, stale_after_sec = lock_stale_after_sec) {
  if (!dir.exists(lock_dir)) {
    return(FALSE)
  }
  lock_age <- timestamp_age_sec(lock_timestamp_from_dir(lock_dir), default = Inf)
  if (!is.finite(lock_age) || lock_age < stale_after_sec) {
    return(FALSE)
  }
  unlink(lock_dir, recursive = TRUE, force = TRUE)
  TRUE
}

with_dir_lock <- function(lock_dir, expr, timeout_sec = 600, poll_sec = 1, stale_after_sec = max(lock_stale_after_sec, timeout_sec + 60L)) {
  expr_sub <- substitute(expr)
  start_time <- Sys.time()
  repeat {
    if (dir.create(lock_dir, showWarnings = FALSE, recursive = FALSE)) {
      write_lock_metadata(lock_dir)
      break
    }
    break_stale_lock_if_needed(lock_dir, stale_after_sec = stale_after_sec)
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    if (is.finite(elapsed) && elapsed >= timeout_sec) {
      stop(sprintf("Timeout while waiting for lock: %s", lock_dir), call. = FALSE)
    }
    Sys.sleep(poll_sec)
  }
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
  eval(expr_sub, envir = parent.frame())
}

assert_file_exists <- function(path, label = basename(path)) {
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", label, path), call. = FALSE)
  }
  invisible(path)
}

ensure_export_path_writable <- function(path) {
  test_file <- file.path(path, ".stage6_write_test.txt")
  ok <- tryCatch({
    writeLines(as.character(Sys.time()), test_file)
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    stop(sprintf("Export path is not writable: %s", path), call. = FALSE)
  }
  unlink(test_file)
}

make_temp_export_path <- function(path) {
  ext <- tools::file_ext(path)
  tempfile(
    pattern = paste0(tools::file_path_sans_ext(basename(path)), "_"),
    tmpdir = dirname(path),
    fileext = if (nzchar(ext)) paste0(".", ext) else ""
  )
}

replace_file_atomically <- function(temp_path, final_path) {
  if (file.exists(final_path)) {
    unlink(final_path)
  }
  renamed <- file.rename(temp_path, final_path)
  if (isTRUE(renamed)) {
    return(TRUE)
  }
  copied <- file.copy(temp_path, final_path, overwrite = TRUE)
  if (isTRUE(copied)) {
    unlink(temp_path)
    return(TRUE)
  }
  FALSE
}

safe_write_csv_atomic <- function(df, path) {
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  readr::write_csv(df, temp_path)
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to place CSV atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_save_rds_atomic <- function(object, path) {
  temp_path <- make_temp_export_path(path)
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  saveRDS(object, temp_path)
  ok <- replace_file_atomically(temp_path, path)
  if (!isTRUE(ok)) {
    stop(sprintf("Failed to place RDS atomically: %s", path), call. = FALSE)
  }
  invisible(path)
}

safe_read_csv <- function(path) {
  safe_read_csv_guess(path)
}

safe_save_plot <- function(plot_object, path, width, height, dpi = 320, bg = "white") {
  tryCatch(
    {
      temp_path <- make_temp_export_path(path)
      on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
      ggplot2::ggsave(
        filename = temp_path,
        plot = plot_object,
        width = width,
        height = height,
        units = "in",
        dpi = dpi,
        bg = bg,
        limitsize = FALSE
      )
      ok <- replace_file_atomically(temp_path, path)
      if (!isTRUE(ok)) {
        stop(sprintf("Failed to place plot atomically: %s", path), call. = FALSE)
      }
      TRUE
    },
    error = function(e) {
      log_step("Failed to save plot: ", path, " :: ", conditionMessage(e), level = "WARN")
      FALSE
    }
  )
}

normalize_runtime_status_df <- function(df) {
  required_cols <- c(
    "timestamp_utc", "worker_id", "stage", "dataset_key", "note", "run_mode",
    "bootstrap_reps_xie", "bootstrap_reps_hsu", "xie_shard_size", "hsu_shard_size", "parallel_cores"
  )
  if (is.null(df) || nrow(df) == 0L) {
    out <- tibble::as_tibble(stats::setNames(rep(list(character()), length(required_cols)), required_cols))
    return(out)
  }
  df <- tibble::as_tibble(df)
  if ("timestamp" %in% names(df) && !"timestamp_utc" %in% names(df)) {
    df$timestamp_utc <- as.character(df$timestamp)
  }
  missing_cols <- setdiff(required_cols, names(df))
  for (nm in missing_cols) {
    df[[nm]] <- NA_character_
  }
  df <- df[, required_cols, drop = FALSE]
  dplyr::mutate(df, dplyr::across(dplyr::everything(), as.character))
}

write_runtime_status <- function(stage, dataset_key = NA_character_, note = NA_character_) {
  status_row <- tibble::tibble(
    timestamp_utc = now_utc_chr(),
    worker_id = as.character(worker_id),
    stage = as.character(stage),
    dataset_key = as.character(dataset_key),
    note = as.character(note),
    run_mode = as.character(run_mode),
    bootstrap_reps_xie = as.character(bootstrap_reps_xie),
    bootstrap_reps_hsu = as.character(bootstrap_reps_hsu),
    xie_shard_size = as.character(xie_shard_size),
    hsu_shard_size = as.character(hsu_shard_size),
    parallel_cores = as.character(parallel_cores)
  )
  with_dir_lock(runtime_lock_dir, {
    existing <- if (file.exists(runtime_status_file)) {
      normalize_runtime_status_df(safe_read_csv_text(runtime_status_file))
    } else {
      normalize_runtime_status_df(NULL)
    }
    combined <- dplyr::bind_rows(existing, status_row)
    safe_write_csv_atomic(combined, runtime_status_file)
  }, timeout_sec = 600, poll_sec = 0.25)
  invisible(TRUE)
}

compact_bind_rows <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (length(x) == 0L) {
    tibble::tibble()
  } else {
    dplyr::bind_rows(x)
  }
}

split_rep_ranges <- function(n_reps, shard_size) {
  if (!is.finite(n_reps) || n_reps <= 0L) {
    return(list())
  }
  idx <- seq_len(n_reps)
  split(idx, ceiling(idx / shard_size))
}

reset_run_state_files <- function() {
  paths_to_remove <- c(
    state_dir,
    observed_dir,
    shards_dir,
    merged_dir,
    shard_registry_file,
    run_manifest_file,
    runtime_status_file,
    progress_log_file
  )
  invisible(lapply(paths_to_remove, function(path) {
    if (file.exists(path) || dir.exists(path)) {
      unlink(path, recursive = TRUE, force = TRUE)
    }
  }))
  dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(observed_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(xie_shards_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(hsu_shards_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(merged_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(heartbeat_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(TRUE)
}

copy_final_outputs_to_canonical <- function(file_names) {
  if (!isTRUE(copy_final_exports_to_canonical)) {
    return(invisible(FALSE))
  }
  if (!nzchar(canonical_export_path) || is.na(canonical_export_path)) {
    return(invisible(FALSE))
  }
  dir.create(canonical_export_path, recursive = TRUE, showWarnings = FALSE)
  for (nm in file_names) {
    src <- file.path(export_path, nm)
    dst <- file.path(canonical_export_path, nm)
    if (file.exists(src)) {
      file.copy(src, dst, overwrite = TRUE)
    }
  }
  invisible(TRUE)
}

make_xie_observed_file <- function(dataset_key) {
  file.path(observed_dir, paste0("xie__", dataset_key, "__observed.rds"))
}

make_hsu_observed_file <- function(dataset_key, family_name) {
  file.path(observed_dir, paste0("hsu__", dataset_key, "__", family_name, "__observed.rds"))
}

make_xie_shard_file <- function(dataset_key, shard_id, n_shards) {
  file.path(xie_shards_dir, sprintf("xie__%s__shard_%03d_of_%03d.rds", dataset_key, shard_id, n_shards))
}

make_hsu_shard_file <- function(dataset_key, family_name, shard_id, n_shards) {
  file.path(hsu_shards_dir, sprintf("hsu__%s__%s__shard_%03d_of_%03d.rds", dataset_key, family_name, shard_id, n_shards))
}

analysis_variant_from_key <- function(dataset_key) {
  dplyr::case_when(
    dataset_key %in% c("PNU", "SNU") ~ "single_cohort",
    dataset_key == "merged__site_free" ~ "site_free",
    dataset_key == "merged__site_adjusted" ~ "site_adjusted",
    TRUE ~ NA_character_
  )
}

source_dataset_from_key <- function(dataset_key) {
  dplyr::case_when(
    dataset_key == "PNU" ~ "PNU",
    dataset_key == "SNU" ~ "SNU",
    dataset_key %in% c("merged__site_free", "merged__site_adjusted") ~ "merged",
    TRUE ~ NA_character_
  )
}

shard_registry_columns <- c(
  "task_id", "task_type", "dataset_key", "source_dataset", "analysis_variant",
  "working_family", "shard_id", "n_shards", "rep_id_start", "rep_id_end", "n_reps",
  "status", "worker_id", "attempt_count", "created_at_utc", "started_at_utc",
  "heartbeat_at_utc", "finished_at_utc", "n_valid_stats", "shard_file", "observed_file", "note"
)

normalize_shard_registry <- function(df) {
  if (is.null(df) || nrow(df) == 0L) {
    out <- tibble::as_tibble(stats::setNames(rep(list(character()), length(shard_registry_columns)), shard_registry_columns))
    return(out)
  }
  df <- tibble::as_tibble(df)
  legacy_map <- c(
    created_at = "created_at_utc",
    started_at = "started_at_utc",
    finished_at = "finished_at_utc",
    heartbeat_at = "heartbeat_at_utc"
  )
  for (old_nm in names(legacy_map)) {
    new_nm <- legacy_map[[old_nm]]
    if (old_nm %in% names(df) && !new_nm %in% names(df)) {
      df[[new_nm]] <- as.character(df[[old_nm]])
    }
  }
  missing_cols <- setdiff(shard_registry_columns, names(df))
  for (nm in missing_cols) {
    df[[nm]] <- NA_character_
  }
  df <- df[, shard_registry_columns, drop = FALSE]
  df <- dplyr::mutate(df, dplyr::across(dplyr::everything(), as.character))
  df$status[is.na(df$status) | !nzchar(df$status)] <- "pending"
  df$attempt_count[is.na(df$attempt_count) | !nzchar(df$attempt_count)] <- "0"
  df
}

read_shard_registry <- function() {
  assert_file_exists(shard_registry_file, "stage6_shard_registry.csv")
  normalize_shard_registry(safe_read_csv_text(shard_registry_file))
}

write_shard_registry <- function(registry_df) {
  safe_write_csv_atomic(normalize_shard_registry(registry_df), shard_registry_file)
  invisible(TRUE)
}

worker_heartbeat_file <- function(id = worker_id) {
  safe_id <- gsub("[^A-Za-z0-9_.-]", "_", as.character(id))
  file.path(heartbeat_dir, paste0(safe_id, ".csv"))
}

touch_worker_heartbeat <- function(task_id = NA_character_, dataset_key = NA_character_, stage = "alive", note = NA_character_) {
  dir.create(heartbeat_dir, recursive = TRUE, showWarnings = FALSE)
  hb <- tibble::tibble(
    timestamp_utc = now_utc_chr(),
    worker_id = as.character(worker_id),
    pid = as.character(Sys.getpid()),
    task_id = as.character(task_id),
    dataset_key = as.character(dataset_key),
    stage = as.character(stage),
    note = as.character(note)
  )
  safe_write_csv_atomic(hb, worker_heartbeat_file(worker_id))
  invisible(TRUE)
}

read_worker_heartbeat <- function(id) {
  path <- worker_heartbeat_file(id)
  if (!file.exists(path)) {
    return(tibble::tibble())
  }
  safe_read_csv_text(path)
}

worker_last_seen_age_sec <- function(id, fallback_timestamp = NA_character_) {
  hb <- read_worker_heartbeat(id)
  hb_ts <- if (nrow(hb) > 0L && "timestamp_utc" %in% names(hb)) hb$timestamp_utc[[1]] else NA_character_
  last_seen <- coalesce_chr_scalar(hb_ts, fallback_timestamp)
  timestamp_age_sec(last_seen, default = Inf)
}

build_stage6_run_manifest <- function() {
  tibble::tibble(
    item = c(
      "run_root_dir", "export_path", "canonical_export_path", "run_mode",
      "bootstrap_reps_xie", "bootstrap_reps_hsu", "xie_shard_size", "hsu_shard_size",
      "candidate_receus_families", "candidate_hsu_families", "hsu_family_set_name",
      "decision_rule_version", "screening_context", "stage6_specification_file",
      "copy_final_exports_to_canonical", "refresh_run_state",
      "task_stale_after_sec", "failed_task_retry_after_sec", "max_task_attempts", "lock_stale_after_sec"
    ),
    value = c(
      run_root_dir, export_path, canonical_export_path, run_mode,
      as.character(bootstrap_reps_xie), as.character(bootstrap_reps_hsu),
      as.character(xie_shard_size), as.character(hsu_shard_size),
      paste(candidate_receus_families, collapse = "|"),
      paste(candidate_hsu_families, collapse = "|"),
      hsu_family_set_name,
      decision_rule_version, screening_context, stage6_specification_file,
      as.character(copy_final_exports_to_canonical), as.character(refresh_run_state),
      as.character(task_stale_after_sec), as.character(failed_task_retry_after_sec),
      as.character(max_task_attempts), as.character(lock_stale_after_sec)
    )
  )
}

# 🔴 Define: Stage 1 ingest wrapped for sharded Stage 6 ===============================
ingest_stage1_assets <- function() {
  required_stage1_files <- c(
    stage1_bundle_file,
    stage1_datasets_file,
    stage1_dataset_registry_file,
    stage1_scaling_registry_file,
    stage1_formula_registry_file,
    stage1_modeling_registry_file,
    stage1_horizon_registry_file,
    stage1_threshold_registry_file,
    stage1_metadata_registry_file
  )
  
  missing_stage1_files <- required_stage1_files[!file.exists(required_stage1_files)]
  if (length(missing_stage1_files) > 0L) {
    stop(
      sprintf(
        "Stage 1 inputs are missing: %s",
        paste(missing_stage1_files, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  ingest_start <- Sys.time()
  log_step("Reading Stage 1 bundles and registries.")
  
  backbone_bundle <- readRDS(stage1_bundle_file)
  analysis_datasets_stage1 <- readRDS(stage1_datasets_file)
  
  dataset_registry <- safe_read_csv(stage1_dataset_registry_file)
  scaling_registry <- safe_read_csv(stage1_scaling_registry_file)
  formula_registry <- safe_read_csv(stage1_formula_registry_file)
  modeling_registry <- safe_read_csv(stage1_modeling_registry_file)
  horizon_registry <- safe_read_csv(stage1_horizon_registry_file)
  threshold_registry <- safe_read_csv(stage1_threshold_registry_file)
  metadata_registry <- safe_read_csv(stage1_metadata_registry_file)
  
  log_step("Stage 1 assets loaded in ", format_elapsed(ingest_start))
  write_runtime_status("stage1_loaded", note = "Stage 1 bundles and registries loaded.")
  
  if (!is.list(analysis_datasets_stage1) || !all(c("PNU", "SNU", "merged") %in% names(analysis_datasets_stage1))) {
    stop("`stage1_analysis_datasets.rds` must contain named datasets: PNU, SNU, merged.", call. = FALSE)
  }
  
  if (!is.list(backbone_bundle) || is.null(backbone_bundle$config) || is.null(backbone_bundle$datasets)) {
    stop("`stage1_backbone_bundle.rds` is malformed or incomplete.", call. = FALSE)
  }
  
  canonical_common_rules <- metadata_registry %>%
    filter(metadata_name == "canonical_common_rules") %>%
    pull(metadata_value)
  
  if (length(canonical_common_rules) == 0L || canonical_common_rules[[1]] != "2.General Model Specifications_🇬🇧ENG.md") {
    stop(
      "Stage 6 requires Stage 1 metadata to identify `2.General Model Specifications_🇬🇧ENG.md` as the canonical governing backbone.",
      call. = FALSE
    )
  }
  
  common_horizons_year <- as.integer(sort(unique(as.integer(unlist(backbone_bundle$config$common_horizons_year)))))
  if (length(common_horizons_year) == 0L && "horizon_year" %in% names(horizon_registry)) {
    common_horizons_year <- as.integer(sort(unique(horizon_registry$horizon_year)))
  }
  if (!identical(common_horizons_year, 1:10)) {
    stop(
      sprintf(
        "Stage 6 requires the inherited horizon grid to be exactly 1:10; found: %s",
        paste(common_horizons_year, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  risk_thresholds <- as.numeric(unlist(backbone_bundle$config$risk_thresholds))
  if (length(risk_thresholds) == 0L && "threshold" %in% names(threshold_registry)) {
    risk_thresholds <- as.numeric(threshold_registry$threshold)
  }
  risk_thresholds <- sort(unique(risk_thresholds))
  if (length(risk_thresholds) == 0L) {
    stop("Stage 6 could not recover the Stage 1 threshold grid.", call. = FALSE)
  }
  
  source_lookup <- backbone_bundle$source_lookup
  if (is.null(source_lookup) || length(source_lookup) == 0L) {
    if (!all(c("dataset", "source_description") %in% names(dataset_registry))) {
      stop("Could not recover Stage 1 source descriptions.", call. = FALSE)
    }
    source_lookup <- as.list(stats::setNames(dataset_registry$source_description, dataset_registry$dataset))
  }
  
  list(
    backbone_bundle = backbone_bundle,
    analysis_datasets_stage1 = analysis_datasets_stage1,
    dataset_registry = dataset_registry,
    scaling_registry = scaling_registry,
    formula_registry = formula_registry,
    modeling_registry = modeling_registry,
    horizon_registry = horizon_registry,
    threshold_registry = threshold_registry,
    metadata_registry = metadata_registry,
    common_horizons_year = common_horizons_year,
    risk_thresholds = risk_thresholds,
    source_lookup = source_lookup
  )
}

dataset_order_full <- c("PNU", "SNU", "merged__site_free", "merged__site_adjusted")
method_order <- c(
  "maller_zhou_dn",
  "maller_zhou_alpha_n",
  "xie_sufficient_followup",
  "receus_weibull",
  "receus_aic",
  "hsu_supscore_weibull_approx",
  "hsu_supscore_lognormal_approx",
  "hsu_supscore_loglogistic_approx",
  hsu_familyset_method_name
)

# 🔴 Define: mathematical and data utilities ===============================
coerce_numeric_text <- function(x) {
  as.numeric(trimws(as.character(x)))
}

check_required_columns <- function(df, required_cols, dataset_name) {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf("[%s] Missing required columns: %s", dataset_name, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
}

safe_log <- function(x) {
  log(pmax(x, .Machine$double.xmin))
}

safe_divide <- function(numerator, denominator) {
  ifelse(is.na(denominator) | denominator == 0, NA_real_, numerator / denominator)
}

clamp_prob <- function(x, lower = 1e-8, upper = 1 - 1e-8) {
  pmin(pmax(x, lower), upper)
}

logit_clip <- function(p) {
  qlogis(clamp_prob(p))
}

compact_character_set <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) {
    NA_character_
  } else {
    paste(sort(x), collapse = "|")
  }
}

append_note <- function(existing_note, addition) {
  ifelse(
    is.na(existing_note) | !nzchar(existing_note),
    addition,
    paste(existing_note, addition, sep = " | ")
  )
}

strip_formula_rhs <- function(formula_text) {
  ifelse(
    is.na(formula_text),
    NA_character_,
    gsub("^\\s*~\\s*", "", as.character(formula_text))
  )
}

extract_list_scalar <- function(x, name, default = NA) {
  if (is.null(x) || is.null(x[[name]]) || length(x[[name]]) == 0L) {
    return(default)
  }
  x[[name]][[1]]
}

extract_frame_cell <- function(df, column, default = NA) {
  if (nrow(df) == 0L || !column %in% names(df)) {
    return(default)
  }
  df[[column]][[1]]
}

safe_exp_param <- function(x, lower = -max_log_param_abs, upper = max_log_param_abs) {
  x <- as.numeric(x)
  if (length(x) == 0L || any(!is.finite(x))) {
    return(rep(NA_real_, length(x)))
  }
  exp(pmin(pmax(x, lower), upper))
}

invalid_density_result <- function(time) {
  list(
    Su = rep(NA_real_, length(time)),
    fu = rep(NA_real_, length(time)),
    params = NULL,
    valid = FALSE
  )
}

run_optim_attempts <- function(start_list, objective_fn, control = list(maxit = 1000, reltol = 1e-10)) {
  start_list <- lapply(start_list, function(x) as.numeric(x))
  candidate_results <- vector("list", length(start_list) * 2L)
  result_index <- 1L
  
  for (start_values in start_list) {
    fit_bfgs <- try(
      stats::optim(par = start_values, fn = objective_fn, method = "BFGS", control = control),
      silent = TRUE
    )
    candidate_results[[result_index]] <- fit_bfgs
    result_index <- result_index + 1L
    
    fit_nm <- try(
      stats::optim(par = start_values, fn = objective_fn, method = "Nelder-Mead", control = control),
      silent = TRUE
    )
    candidate_results[[result_index]] <- fit_nm
    result_index <- result_index + 1L
  }
  
  valid_results <- candidate_results[!vapply(candidate_results, inherits, logical(1), what = "try-error")]
  valid_results <- valid_results[vapply(valid_results, function(x) is.list(x) && is.finite(x$value), logical(1))]
  
  if (length(valid_results) == 0L) {
    return(NULL)
  }
  
  valid_results[[which.min(vapply(valid_results, function(x) x$value, numeric(1)))]]
}

run_bootstrap_chunks <- function(n_reps, seed, label, worker_fun) {
  if (n_reps <= 0L) {
    return(numeric(0))
  }
  
  results <- rep(NA_real_, n_reps)
  chunk_ids <- split(seq_len(n_reps), ceiling(seq_len(n_reps) / bootstrap_chunk_size))
  n_chunks <- length(chunk_ids)
  use_parallel <- isTRUE(allow_parallel_bootstrap) && .Platform$OS.type != "windows" && isTRUE(parallel_cores > 1L)
  
  for (chunk_index in seq_along(chunk_ids)) {
    idx <- chunk_ids[[chunk_index]]
    log_step(
      label,
      ": bootstrap chunk ",
      chunk_index,
      "/",
      n_chunks,
      " (",
      min(idx),
      "-",
      max(idx),
      ")"
    )
    
    run_one <- function(rep_id) {
      set.seed(seed + rep_id)
      tryCatch(as.numeric(worker_fun(rep_id)), error = function(e) NA_real_)
    }
    
    chunk_result <- if (use_parallel && length(idx) > 1L) {
      unlist(parallel::mclapply(idx, run_one, mc.cores = min(parallel_cores, length(idx))), use.names = FALSE)
    } else {
      vapply(idx, run_one, numeric(1))
    }
    
    results[idx] <- chunk_result
    log_step(
      label,
      ": completed chunk ",
      chunk_index,
      "/",
      n_chunks,
      " with ",
      sum(is.finite(chunk_result)),
      "/",
      length(chunk_result),
      " valid statistics"
    )
    write_runtime_status("bootstrap_running", note = paste(label, "chunk", chunk_index, "of", n_chunks))
    invisible(gc(verbose = FALSE))
  }
  
  results
}

get_scaling_row <- function(scaling_registry, dataset_name) {
  scaling_row <- scaling_registry[
    scaling_registry$dataset == dataset_name &
      scaling_registry$variable == "age_exact_entry" &
      scaling_registry$scaled_variable == "age_s",
    ,
    drop = FALSE
  ]
  
  if (nrow(scaling_row) == 0L) {
    scaling_row <- scaling_registry[scaling_registry$dataset == dataset_name, , drop = FALSE]
  }
  
  if (nrow(scaling_row) == 0L) {
    NULL
  } else {
    scaling_row[1, , drop = FALSE]
  }
}

select_formula_rhs <- function(formula_registry, dataset_name, site_branch = NULL, formula_name = NULL) {
  out <- formula_registry[formula_registry$dataset == dataset_name, , drop = FALSE]
  
  if (!is.null(site_branch) && "site_branch" %in% names(out)) {
    out <- out[out$site_branch == site_branch, , drop = FALSE]
  }
  
  if (!is.null(formula_name) && "formula_name" %in% names(out)) {
    out <- out[out$formula_name == formula_name, , drop = FALSE]
  }
  
  unique(as.character(out$formula_rhs))
}

limit_formula_set <- function(x, dataset_key) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  
  if (length(x) == 0L) {
    return(x)
  }
  
  if (is.finite(max_hsu_formula_count) && length(x) > max_hsu_formula_count) {
    log_step(
      "HSU incidence formula set for ",
      dataset_key,
      " truncated from ",
      length(x),
      " to ",
      max_hsu_formula_count,
      " because run_mode=",
      run_mode
    )
    x <- x[seq_len(max_hsu_formula_count)]
  }
  
  x
}

hsu_method_name <- function(family) {
  switch(
    as.character(family),
    weibull = "hsu_supscore_weibull_approx",
    lognormal = "hsu_supscore_lognormal_approx",
    loglogistic = "hsu_supscore_loglogistic_approx",
    stop(sprintf("Unsupported HSU family for method naming: %s", family), call. = FALSE)
  )
}

# 🔴 Build: dataset normalization helpers ===============================
prepare_stage6_dataset <- function(df, dataset_name, scaling_registry) {
  required_cols <- c("id", "site", "sex_num", "age_exact_entry", "days_followup", "status_num")
  check_required_columns(df, required_cols, dataset_name)
  
  df <- df %>%
    mutate(
      id = trimws(as.character(id)),
      site = trimws(as.character(site)),
      sex_num = as.integer(coerce_numeric_text(sex_num)),
      age_exact_entry = coerce_numeric_text(age_exact_entry),
      days_followup = coerce_numeric_text(days_followup),
      status_num = as.integer(coerce_numeric_text(status_num))
    )
  
  if (nrow(df) == 0L) {
    stop(sprintf("[%s] Dataset has zero rows.", dataset_name), call. = FALSE)
  }
  
  if (anyNA(df[required_cols])) {
    stop(sprintf("[%s] Missing values detected in Stage 6 backbone columns.", dataset_name), call. = FALSE)
  }
  
  if (any(df$id == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `id` values detected.", dataset_name), call. = FALSE)
  }
  
  if (any(df$site == "", na.rm = TRUE)) {
    stop(sprintf("[%s] Blank `site` values detected.", dataset_name), call. = FALSE)
  }
  
  if (any(!df$sex_num %in% c(0L, 1L))) {
    stop(sprintf("[%s] `sex_num` must be coded as 0/1.", dataset_name), call. = FALSE)
  }
  
  if (any(!df$status_num %in% c(0L, 1L, 2L))) {
    stop(sprintf("[%s] `status_num` must be coded as 0/1/2.", dataset_name), call. = FALSE)
  }
  
  if (any(df$days_followup < 0)) {
    stop(sprintf("[%s] Negative `days_followup` values detected.", dataset_name), call. = FALSE)
  }
  
  if (!"unique_person_id" %in% names(df)) {
    df$unique_person_id <- paste(df$site, df$id, sep = "_")
  } else {
    df$unique_person_id <- trimws(as.character(df$unique_person_id))
  }
  
  if (nrow(df) != dplyr::n_distinct(df$unique_person_id)) {
    stop(sprintf("[%s] `site + id` is not unique within dataset.", dataset_name), call. = FALSE)
  }
  
  if (!"time_year" %in% names(df)) {
    df$time_year <- df$days_followup / 365.25
  } else {
    df$time_year <- coerce_numeric_text(df$time_year)
  }
  
  if (!"event_main" %in% names(df)) {
    df$event_main <- as.integer(df$status_num == 1L)
  } else {
    df$event_main <- as.integer(coerce_numeric_text(df$event_main))
  }
  
  if (!"right_censor_flag" %in% names(df)) {
    df$right_censor_flag <- as.integer(df$status_num == 0L)
  } else {
    df$right_censor_flag <- as.integer(coerce_numeric_text(df$right_censor_flag))
  }
  
  if (!"remission_flag" %in% names(df)) {
    df$remission_flag <- as.integer(df$status_num == 2L)
  } else {
    df$remission_flag <- as.integer(coerce_numeric_text(df$remission_flag))
  }
  
  if (!"censor_main" %in% names(df)) {
    df$censor_main <- as.integer(df$status_num %in% c(0L, 2L))
  } else {
    df$censor_main <- as.integer(coerce_numeric_text(df$censor_main))
  }
  
  if (!"sex_label" %in% names(df)) {
    df$sex_label <- factor(ifelse(df$sex_num == 0L, "Female", "Male"), levels = c("Female", "Male"))
  } else {
    df$sex_label <- factor(as.character(df$sex_label), levels = c("Female", "Male"))
  }
  
  if (!"status_label" %in% names(df)) {
    df$status_label <- factor(
      dplyr::case_when(
        df$status_num == 0L ~ "right_censoring",
        df$status_num == 2L ~ "remission",
        TRUE ~ "transition"
      ),
      levels = c("right_censoring", "remission", "transition")
    )
  } else {
    df$status_label <- factor(as.character(df$status_label), levels = c("right_censoring", "remission", "transition"))
  }
  
  if (!"age_s" %in% names(df) || anyNA(df$age_s)) {
    scaling_row <- get_scaling_row(scaling_registry, dataset_name)
    
    if (!is.null(scaling_row) &&
        all(c("center_mean", "scale_two_sd") %in% names(scaling_row)) &&
        is.finite(scaling_row$center_mean[[1]]) &&
        is.finite(scaling_row$scale_two_sd[[1]]) &&
        scaling_row$scale_two_sd[[1]] > 0) {
      df$age_s <- (df$age_exact_entry - scaling_row$center_mean[[1]]) / scaling_row$scale_two_sd[[1]]
    } else {
      age_mean <- mean(df$age_exact_entry)
      age_sd <- stats::sd(df$age_exact_entry)
      if (is.na(age_sd) || age_sd <= 0) {
        stop(sprintf("[%s] `age_exact_entry` must have positive standard deviation.", dataset_name), call. = FALSE)
      }
      df$age_s <- (df$age_exact_entry - age_mean) / (2 * age_sd)
    }
  } else {
    df$age_s <- coerce_numeric_text(df$age_s)
  }
  
  if (dataset_name %in% c("PNU", "SNU") && dplyr::n_distinct(df$site) != 1L) {
    stop(sprintf("[%s] Single-cohort Stage 6 input must contain exactly one site.", dataset_name), call. = FALSE)
  }
  
  if (dataset_name == "merged" && dplyr::n_distinct(df$site) < 2L) {
    stop("[merged] Merged Stage 6 input must contain at least two site levels.", call. = FALSE)
  }
  
  df
}

# 🔴 Build: Kaplan-Meier evaluators ===============================
make_km_fit <- function(time, event) {
  survival::survfit(survival::Surv(time, event) ~ 1)
}

km_survival_from_fit <- function(km_fit, times) {
  times <- pmax(as.numeric(times), 0)
  if (length(times) == 0L) {
    return(numeric(0))
  }
  if (length(km_fit$time) == 0L) {
    return(rep(1, length(times)))
  }
  position <- findInterval(times, km_fit$time)
  surv_out <- rep(1, length(times))
  positive_index <- position > 0L
  surv_out[positive_index] <- km_fit$surv[position[positive_index]]
  as.numeric(surv_out)
}

km_survival_at <- function(time, event, times, km_fit = NULL) {
  if (is.null(km_fit)) {
    km_fit <- make_km_fit(time, event)
  }
  km_survival_from_fit(km_fit, times)
}

km_cdf_at <- function(time, event, times, km_fit = NULL) {
  1 - km_survival_at(time, event, times, km_fit)
}

compute_mz_extreme_summaries <- function(time, event) {
  time <- as.numeric(time)
  event <- as.integer(event)
  n <- length(time)
  t_n <- max(time)
  
  if (!any(event == 1L)) {
    return(list(t_n = t_n, t_k = NA_real_, lower_bound = NA_real_, N_n = NA_integer_, q_n = NA_real_, alpha_n = NA_real_))
  }
  
  t_k <- max(time[event == 1L])
  lower_bound <- 2 * t_k - t_n
  N_n <- sum(event == 1L & time > lower_bound & time <= t_k)
  q_n <- safe_divide(N_n, n)
  alpha_n <- if (is.na(q_n)) NA_real_ else (1 - q_n)^n
  
  list(
    t_n = t_n,
    t_k = t_k,
    lower_bound = lower_bound,
    N_n = as.integer(N_n),
    q_n = as.numeric(q_n),
    alpha_n = as.numeric(alpha_n)
  )
}

compute_xie_statistic <- function(time, event) {
  time <- as.numeric(time)
  event <- as.integer(event)
  time <- pmax(time, 0)
  km_fit <- make_km_fit(time, event)
  
  t_n <- max(time)
  p_hat_n <- as.numeric(km_cdf_at(time, event, t_n, km_fit))
  
  if (!any(event == 1L)) {
    return(list(t_n = t_n, t_k = NA_real_, epsilon = NA_real_, p_hat_n = p_hat_n, p_hat_G = NA_real_, T_n = NA_real_, degenerate = TRUE))
  }
  
  t_k <- max(time[event == 1L])
  epsilon <- if (2 * (t_n - t_k) < t_n) {
    (9 / 8) * t_n - 0.25 * t_k
  } else {
    t_n
  }
  
  if (2 * t_k < t_n) {
    return(list(t_n = t_n, t_k = t_k, epsilon = epsilon, p_hat_n = p_hat_n, p_hat_G = p_hat_n, T_n = 0, degenerate = TRUE))
  }
  
  t_1 <- max(0, t_n - epsilon)
  t_2 <- max(0, t_n - epsilon / 2)
  f_t1 <- as.numeric(km_cdf_at(time, event, t_1, km_fit))
  f_t2 <- as.numeric(km_cdf_at(time, event, t_2, km_fit))
  denominator <- 2 * f_t2 - f_t1 - p_hat_n
  
  p_hat_G <- if (!is.finite(denominator) || denominator <= 0) {
    p_hat_n
  } else {
    f_t1 + (f_t2 - f_t1)^2 / denominator
  }
  
  if (!is.finite(p_hat_G) || p_hat_G < p_hat_n) {
    p_hat_G <- p_hat_n
  }
  p_hat_G <- min(max(p_hat_G, p_hat_n), 1)
  
  list(
    t_n = t_n,
    t_k = t_k,
    epsilon = epsilon,
    p_hat_n = p_hat_n,
    p_hat_G = p_hat_G,
    T_n = max(0, p_hat_G - p_hat_n),
    degenerate = FALSE
  )
}

run_xie_bootstrap <- function(df, bootstrap_reps, seed, label) {
  time <- as.numeric(df$time_year)
  event <- as.integer(df$event_main)
  observed <- compute_xie_statistic(time, event)
  
  if (is.na(observed$T_n)) {
    return(list(
      observed = observed,
      p_centered = NA_real_,
      p_raw = NA_real_,
      critical_centered = NA_real_,
      bootstrap_stats = rep(NA_real_, bootstrap_reps)
    ))
  }
  
  if (isTRUE(observed$degenerate) && observed$T_n == 0) {
    return(list(
      observed = observed,
      p_centered = 1,
      p_raw = 1,
      critical_centered = 0,
      bootstrap_stats = rep(0, bootstrap_reps)
    ))
  }
  
  n <- nrow(df)
  bootstrap_stats <- run_bootstrap_chunks(
    n_reps = bootstrap_reps,
    seed = seed,
    label = label,
    worker_fun = function(rep_id) {
      draw_index <- sample.int(n = n, size = n, replace = TRUE)
      compute_xie_statistic(time[draw_index], event[draw_index])$T_n
    }
  )
  
  centered_stats <- bootstrap_stats - observed$T_n
  n_valid_centered <- sum(is.finite(centered_stats))
  n_valid_raw <- sum(is.finite(bootstrap_stats))
  
  critical_centered <- if (n_valid_centered > 0L) {
    as.numeric(stats::quantile(centered_stats, probs = 1 - alpha_screening, na.rm = TRUE, type = 8, names = FALSE))
  } else {
    NA_real_
  }
  
  p_centered <- if (n_valid_centered > 0L) {
    (sum(centered_stats >= observed$T_n, na.rm = TRUE) + 1) / (n_valid_centered + 1)
  } else {
    NA_real_
  }
  
  p_raw <- if (n_valid_raw > 0L) {
    (sum(bootstrap_stats >= observed$T_n, na.rm = TRUE) + 1) / (n_valid_raw + 1)
  } else {
    NA_real_
  }
  
  list(
    observed = observed,
    p_centered = p_centered,
    p_raw = p_raw,
    critical_centered = critical_centered,
    bootstrap_stats = bootstrap_stats
  )
}

# 🔴 Build: intercept-only cure likelihoods for RECeUS ===============================
get_survreg_start <- function(time, event, dist_name) {
  df_tmp <- data.frame(time_survreg = pmax(time, tiny_time), event_survreg = event)
  fit <- try(
    suppressWarnings(
      survival::survreg(
        survival::Surv(time_survreg, event_survreg) ~ 1,
        data = df_tmp,
        dist = dist_name
      )
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) {
    return(NULL)
  }
  
  list(intercept = unname(stats::coef(fit)[1]), survreg_scale = fit$scale)
}

family_surv_density <- function(time, family, par_positive) {
  time <- pmax(as.numeric(time), tiny_time)
  par_positive <- as.numeric(par_positive)
  
  if (length(par_positive) == 0L || any(!is.finite(par_positive))) {
    return(invalid_density_result(time))
  }
  
  if (family == "exponential") {
    rate <- safe_exp_param(par_positive[1])[[1]]
    if (!is.finite(rate) || rate <= 0) {
      return(invalid_density_result(time))
    }
    Su <- exp(-rate * time)
    fu <- rate * Su
    params <- list(rate = rate)
  } else if (family == "weibull") {
    shape <- safe_exp_param(par_positive[1])[[1]]
    scale <- safe_exp_param(par_positive[2])[[1]]
    if (!is.finite(shape) || !is.finite(scale) || shape <= 0 || scale <= 0) {
      return(invalid_density_result(time))
    }
    z <- (time / scale)^shape
    Su <- exp(-z)
    fu <- (shape / scale) * (time / scale)^(shape - 1) * Su
    params <- list(shape = shape, scale = scale)
  } else if (family == "gamma") {
    shape <- safe_exp_param(par_positive[1])[[1]]
    rate <- safe_exp_param(par_positive[2])[[1]]
    if (!is.finite(shape) || !is.finite(rate) || shape <= 0 || rate <= 0) {
      return(invalid_density_result(time))
    }
    gamma_eval <- try(
      suppressWarnings(
        list(
          Su = pgamma(time, shape = shape, rate = rate, lower.tail = FALSE),
          fu = dgamma(time, shape = shape, rate = rate)
        )
      ),
      silent = TRUE
    )
    if (inherits(gamma_eval, "try-error")) {
      return(invalid_density_result(time))
    }
    Su <- gamma_eval$Su
    fu <- gamma_eval$fu
    params <- list(shape = shape, rate = rate)
  } else if (family == "loglogistic") {
    shape <- safe_exp_param(par_positive[1])[[1]]
    scale <- safe_exp_param(par_positive[2])[[1]]
    if (!is.finite(shape) || !is.finite(scale) || shape <= 0 || scale <= 0) {
      return(invalid_density_result(time))
    }
    z <- (time / scale)^shape
    Su <- 1 / (1 + z)
    fu <- (shape / scale) * (time / scale)^(shape - 1) / (1 + z)^2
    params <- list(shape = shape, scale = scale)
  } else {
    stop(sprintf("Unsupported family: %s", family), call. = FALSE)
  }
  
  if (any(!is.finite(Su)) || any(!is.finite(fu)) || any(Su < 0) || any(Su > 1) || any(fu < 0)) {
    return(invalid_density_result(time))
  }
  
  list(
    Su = pmin(pmax(Su, .Machine$double.xmin), 1),
    fu = pmax(fu, .Machine$double.xmin),
    params = params,
    valid = TRUE
  )
}

make_intercept_starts <- function(time, event, family, cure_start) {
  event_time <- time[event == 1L]
  mean_time <- if (length(event_time) > 0L) mean(event_time) else mean(time)
  median_time <- if (length(event_time) > 0L) stats::median(event_time) else stats::median(time)
  
  mean_time <- max(mean_time, tiny_time)
  median_time <- max(median_time, tiny_time)
  
  cure_candidates <- unique(c(
    logit_clip(cure_start),
    logit_clip(0.05),
    logit_clip(0.10),
    logit_clip(min(0.75, max(0.15, cure_start * 1.25))),
    logit_clip(max(0.01, cure_start / 2))
  ))
  
  if (family == "exponential") {
    rate0 <- if (sum(event) > 0L) sum(event) / sum(pmax(time, tiny_time)) else 1 / mean_time
    parameter_candidates <- list(c(log(rate0)), c(log(1 / median_time)), c(log(1 / mean_time)))
  } else if (family == "weibull") {
    sr <- get_survreg_start(time, event, "weibull")
    shape0 <- if (!is.null(sr)) 1 / sr$survreg_scale else 1
    scale0 <- if (!is.null(sr)) exp(sr$intercept) else median_time
    parameter_candidates <- list(c(log(shape0), log(scale0)), c(log(1), log(scale0)), c(log(1.5), log(scale0)))
  } else if (family == "gamma") {
    shape0 <- 1
    rate0 <- shape0 / mean_time
    parameter_candidates <- list(c(log(shape0), log(rate0)), c(log(0.5), log(0.5 / mean_time)), c(log(2), log(2 / mean_time)))
  } else if (family == "loglogistic") {
    sr <- get_survreg_start(time, event, "loglogistic")
    shape0 <- if (!is.null(sr)) 1 / sr$survreg_scale else 1
    scale0 <- if (!is.null(sr)) exp(sr$intercept) else median_time
    parameter_candidates <- list(c(log(shape0), log(scale0)), c(log(1), log(scale0)), c(log(1.5), log(scale0)))
  } else {
    stop(sprintf("Unsupported family: %s", family), call. = FALSE)
  }
  
  start_list <- vector("list", length(cure_candidates) * length(parameter_candidates))
  idx <- 1L
  for (cure_eta in cure_candidates) {
    for (par_candidate in parameter_candidates) {
      start_list[[idx]] <- c(cure_eta, par_candidate)
      idx <- idx + 1L
    }
  }
  start_list
}

negloglik_mix_intercept <- function(par, time, event, family) {
  if (any(!is.finite(par))) {
    return(1e20)
  }
  
  cure_fraction <- clamp_prob(plogis(par[1]))
  family_eval <- family_surv_density(time, family, par[-1])
  
  if (!isTRUE(family_eval$valid)) {
    return(1e20)
  }
  
  Su <- family_eval$Su
  fu <- family_eval$fu
  
  loglik_vector <- ifelse(
    event == 1L,
    safe_log(1 - cure_fraction) + safe_log(fu),
    safe_log(cure_fraction + (1 - cure_fraction) * Su)
  )
  
  if (any(!is.finite(loglik_vector))) {
    return(1e20)
  }
  
  -sum(loglik_vector)
}

fit_mixture_intercept_family <- function(time, event, family) {
  time <- pmax(as.numeric(time), tiny_time)
  event <- as.integer(event)
  tau <- max(time)
  
  km_fit <- make_km_fit(time, event)
  cure_start <- clamp_prob(km_survival_at(time, event, tau, km_fit), lower = 0.001, upper = 0.95)
  start_list <- make_intercept_starts(time, event, family, cure_start)
  objective_fn <- function(par) negloglik_mix_intercept(par, time, event, family)
  
  opt <- run_optim_attempts(
    start_list = start_list,
    objective_fn = objective_fn,
    control = list(maxit = optim_maxit, reltol = optim_reltol)
  )
  
  if (is.null(opt)) {
    return(list(
      converged = FALSE,
      family = family,
      loglik = NA_real_,
      aic = NA_real_,
      cure_fraction_hat = NA_real_,
      susceptible_survival_tau = NA_real_,
      overall_survival_tau = NA_real_,
      receus_ratio_hat = NA_real_,
      tau_year = tau,
      positive_params = NULL,
      optim_par = NA_real_,
      optim_convergence = NA_integer_,
      note = "No valid optimization result."
    ))
  }
  
  cure_fraction_hat <- clamp_prob(plogis(opt$par[1]))
  tau_eval <- family_surv_density(tau, family, opt$par[-1])
  
  if (!isTRUE(tau_eval$valid)) {
    return(list(
      converged = FALSE,
      family = family,
      loglik = NA_real_,
      aic = NA_real_,
      cure_fraction_hat = NA_real_,
      susceptible_survival_tau = NA_real_,
      overall_survival_tau = NA_real_,
      receus_ratio_hat = NA_real_,
      tau_year = tau,
      positive_params = NULL,
      optim_par = opt$par,
      optim_convergence = opt$convergence,
      note = "Optimization finished but tau survival evaluation was invalid."
    ))
  }
  
  susceptible_survival_tau <- as.numeric(tau_eval$Su[1])
  overall_survival_tau <- cure_fraction_hat + (1 - cure_fraction_hat) * susceptible_survival_tau
  receus_ratio_hat <- safe_divide(susceptible_survival_tau, overall_survival_tau)
  if (is.finite(receus_ratio_hat)) {
    receus_ratio_hat <- min(max(receus_ratio_hat, 0), 1)
  }
  
  loglik <- -opt$value
  npar <- length(opt$par)
  
  list(
    converged = isTRUE(opt$convergence == 0L) && is.finite(loglik),
    family = family,
    loglik = loglik,
    aic = 2 * npar - 2 * loglik,
    cure_fraction_hat = cure_fraction_hat,
    susceptible_survival_tau = susceptible_survival_tau,
    overall_survival_tau = overall_survival_tau,
    receus_ratio_hat = receus_ratio_hat,
    tau_year = tau,
    positive_params = tau_eval$params,
    optim_par = opt$par,
    optim_convergence = opt$convergence,
    note = NA_character_
  )
}

fit_no_cure_exponential <- function(time, event) {
  time <- pmax(as.numeric(time), tiny_time)
  event <- as.integer(event)
  
  rate_hat <- if (sum(event) > 0L) {
    sum(event) / sum(time)
  } else {
    1 / max(mean(time), tiny_time)
  }
  
  loglik <- sum(event * log(rate_hat) - rate_hat * time)
  list(rate_hat = rate_hat, loglik = loglik, converged = is.finite(rate_hat) && is.finite(loglik))
}

run_maller_zhou_presence <- function(df) {
  time <- as.numeric(df$time_year)
  event <- as.integer(df$event_main)
  km_fit <- make_km_fit(time, event)
  tau <- max(time)
  km_tail_survival <- as.numeric(km_survival_at(time, event, tau, km_fit))
  extreme_summary <- compute_mz_extreme_summaries(time, event)
  
  if (sum(event) == 0L) {
    return(list(
      result = list(dn = 0, p_value = 1, cure_fraction_hat = 1, km_tail_survival = km_tail_survival, alpha_n = extreme_summary$alpha_n, q_n = extreme_summary$q_n),
      fit_object = list(null = NULL, alt = NULL, extreme_summary = extreme_summary)
    ))
  }
  
  null_fit <- fit_no_cure_exponential(time, event)
  alt_fit <- fit_mixture_intercept_family(time, event, "exponential")
  
  dn <- if (isTRUE(alt_fit$converged) && isTRUE(null_fit$converged)) {
    max(0, 2 * (alt_fit$loglik - null_fit$loglik))
  } else {
    NA_real_
  }
  
  p_value <- if (is.na(dn)) {
    NA_real_
  } else if (dn > 0) {
    0.5 * stats::pchisq(dn, df = 1, lower.tail = FALSE)
  } else {
    1
  }
  
  list(
    result = list(
      dn = dn,
      p_value = p_value,
      cure_fraction_hat = alt_fit$cure_fraction_hat,
      km_tail_survival = km_tail_survival,
      alpha_n = extreme_summary$alpha_n,
      q_n = extreme_summary$q_n
    ),
    fit_object = list(null = null_fit, alt = alt_fit, extreme_summary = extreme_summary)
  )
}

run_receus_suite <- function(df) {
  time <- as.numeric(df$time_year)
  event <- as.integer(df$event_main)
  
  fit_list <- stats::setNames(
    lapply(candidate_receus_families, function(family_name) fit_mixture_intercept_family(time, event, family_name)),
    candidate_receus_families
  )
  
  candidate_table <- bind_rows(lapply(names(fit_list), function(family_name) {
    fit_obj <- fit_list[[family_name]]
    
    tibble::tibble(
      family = family_name,
      converged = isTRUE(fit_obj$converged),
      loglik = extract_list_scalar(fit_obj, "loglik", NA_real_),
      aic = extract_list_scalar(fit_obj, "aic", NA_real_),
      cure_fraction_hat = extract_list_scalar(fit_obj, "cure_fraction_hat", NA_real_),
      susceptible_survival_tau = extract_list_scalar(fit_obj, "susceptible_survival_tau", NA_real_),
      overall_survival_tau = extract_list_scalar(fit_obj, "overall_survival_tau", NA_real_),
      receus_ratio_hat = extract_list_scalar(fit_obj, "receus_ratio_hat", NA_real_),
      tau_year = extract_list_scalar(fit_obj, "tau_year", NA_real_),
      receus_flag = NA_character_,
      selected_by_aic = FALSE,
      note = extract_list_scalar(fit_obj, "note", NA_character_)
    )
  }))
  
  candidate_table <- candidate_table %>%
    mutate(
      receus_flag = vapply(
        seq_len(n()),
        function(i) classify_receus_flag(cure_fraction_hat = cure_fraction_hat[[i]], receus_ratio_hat = receus_ratio_hat[[i]]),
        character(1)
      )
    )
  
  converged_candidates <- candidate_table %>% filter(converged, is.finite(aic))
  best_family <- if (nrow(converged_candidates) > 0L) converged_candidates$family[[which.min(converged_candidates$aic)]] else NA_character_
  
  candidate_table <- candidate_table %>%
    mutate(selected_by_aic = family == best_family)
  
  list(
    fits = fit_list,
    candidate_table = candidate_table,
    best_family = best_family,
    best_fit = if (!is.na(best_family)) fit_list[[best_family]] else NULL,
    weibull_fit = fit_list[["weibull"]]
  )
}

# 🔴 Build: multi-family HSU pragmatic approximation routines ===============================
survreg_dist_lookup <- c(
  weibull = "weibull",
  lognormal = "lognormal",
  loglogistic = "loglogistic"
)

aft_surv_density <- function(time, family, mu, log_aux) {
  time <- pmax(as.numeric(time), tiny_time)
  mu <- as.numeric(mu)
  if (length(mu) != length(time)) {
    stop("`mu` must have the same length as `time`.", call. = FALSE)
  }
  
  aux <- safe_exp_param(log_aux)[[1]]
  if (!is.finite(aux) || aux <= 0) {
    return(invalid_density_result(time))
  }
  
  log_time <- safe_log(time)
  
  if (family == "weibull") {
    z <- (log_time - mu) / aux
    w <- exp(pmin(z, max_log_param_abs))
    Su <- exp(-w)
    fu <- (w * Su) / (aux * time)
    params <- list(sigma = aux)
  } else if (family == "lognormal") {
    z <- (log_time - mu) / aux
    Su <- stats::pnorm(z, lower.tail = FALSE)
    fu <- stats::dnorm(z) / (aux * time)
    params <- list(sdlog = aux)
  } else if (family == "loglogistic") {
    z <- (log_time - mu) / aux
    Su <- stats::plogis(z, lower.tail = FALSE)
    fu <- stats::dlogis(z) / (aux * time)
    params <- list(sigma = aux)
  } else {
    stop(sprintf("Unsupported AFT working family: %s", family), call. = FALSE)
  }
  
  if (any(!is.finite(Su)) || any(!is.finite(fu)) || any(Su < 0) || any(Su > 1) || any(fu < 0)) {
    return(invalid_density_result(time))
  }
  
  list(
    Su = pmin(pmax(Su, .Machine$double.xmin), 1),
    fu = pmax(fu, .Machine$double.xmin),
    params = params,
    valid = TRUE
  )
}

negloglik_no_cure_aft <- function(par, X, time, event, family) {
  p <- ncol(X)
  beta <- par[seq_len(p)]
  log_aux <- par[p + 1L]
  
  if (any(!is.finite(beta)) || !is.finite(log_aux) || abs(log_aux) > max_log_param_abs) {
    return(1e20)
  }
  
  mu <- drop(X %*% beta)
  if (any(!is.finite(mu))) {
    return(1e20)
  }
  
  eval_family <- aft_surv_density(time = time, family = family, mu = mu, log_aux = log_aux)
  if (!isTRUE(eval_family$valid)) {
    return(1e20)
  }
  
  loglik_vector <- ifelse(
    event == 1L,
    safe_log(eval_family$fu),
    safe_log(eval_family$Su)
  )
  
  if (any(!is.finite(loglik_vector))) {
    return(1e20)
  }
  
  -sum(loglik_vector)
}

fit_no_cure_aft <- function(df, latency_formula, family) {
  df_model <- df
  if ("site" %in% names(df_model)) {
    df_model$site <- factor(df_model$site)
  }
  
  time <- pmax(as.numeric(df_model$time_year), tiny_time)
  event <- as.integer(df_model$event_main)
  X <- model.matrix(stats::as.formula(latency_formula), data = df_model)
  
  beta_start <- rep(0, ncol(X))
  names(beta_start) <- colnames(X)
  log_aux_start <- 0
  
  df_survreg <- df_model
  df_survreg$time_survreg <- time
  df_survreg$event_survreg <- event
  
  survreg_fit <- try(
    suppressWarnings(
      survival::survreg(
        stats::as.formula(paste("survival::Surv(time_survreg, event_survreg)", latency_formula)),
        data = df_survreg,
        dist = unname(survreg_dist_lookup[[family]])
      )
    ),
    silent = TRUE
  )
  
  if (!inherits(survreg_fit, "try-error")) {
    survreg_beta <- stats::coef(survreg_fit)
    beta_start[names(survreg_beta)] <- survreg_beta
    log_aux_start <- log(survreg_fit$scale)
  }
  
  start_main <- c(beta_start, log_aux_start)
  start_zero <- c(rep(0, ncol(X)), 0)
  
  objective_fn <- function(par) negloglik_no_cure_aft(par, X, time, event, family)
  
  opt <- run_optim_attempts(
    start_list = list(start_main, start_zero),
    objective_fn = objective_fn,
    control = list(maxit = optim_maxit, reltol = optim_reltol)
  )
  
  if (is.null(opt)) {
    return(list(
      converged = FALSE,
      latency_formula = latency_formula,
      family = family,
      beta_hat = NA_real_,
      aux_hat = NA_real_,
      log_aux_hat = NA_real_,
      loglik = NA_real_,
      aic = NA_real_,
      optim_par = NA_real_,
      optim_convergence = NA_integer_
    ))
  }
  
  beta_hat <- opt$par[seq_len(ncol(X))]
  names(beta_hat) <- colnames(X)
  log_aux_hat <- opt$par[ncol(X) + 1L]
  aux_hat <- safe_exp_param(log_aux_hat)[[1]]
  loglik <- -opt$value
  
  list(
    converged = isTRUE(opt$convergence == 0L) && is.finite(loglik) && is.finite(aux_hat),
    latency_formula = latency_formula,
    family = family,
    beta_hat = beta_hat,
    aux_hat = aux_hat,
    log_aux_hat = log_aux_hat,
    loglik = loglik,
    aic = 2 * length(opt$par) - 2 * loglik,
    optim_par = opt$par,
    optim_convergence = opt$convergence
  )
}

negloglik_mix_aft <- function(par, Z, X, time, event, family) {
  q <- ncol(Z)
  p <- ncol(X)
  
  gamma <- par[seq_len(q)]
  beta <- par[q + seq_len(p)]
  log_aux <- par[q + p + 1L]
  
  if (any(!is.finite(gamma)) || any(!is.finite(beta)) || !is.finite(log_aux) || abs(log_aux) > max_log_param_abs) {
    return(1e20)
  }
  
  eta <- drop(Z %*% gamma)
  mu <- drop(X %*% beta)
  if (any(!is.finite(eta)) || any(!is.finite(mu))) {
    return(1e20)
  }
  
  cure_i <- clamp_prob(plogis(eta))
  eval_family <- aft_surv_density(time = time, family = family, mu = mu, log_aux = log_aux)
  if (!isTRUE(eval_family$valid)) {
    return(1e20)
  }
  
  loglik_vector <- ifelse(
    event == 1L,
    safe_log(1 - cure_i) + safe_log(eval_family$fu),
    safe_log(cure_i + (1 - cure_i) * eval_family$Su)
  )
  
  if (any(!is.finite(loglik_vector))) {
    return(1e20)
  }
  
  -sum(loglik_vector)
}

fit_mixture_cure_aft <- function(df, incidence_formula, latency_formula, family, null_beta_start = NULL, null_log_aux_start = NULL) {
  df_model <- df
  if ("site" %in% names(df_model)) {
    df_model$site <- factor(df_model$site)
  }
  
  time <- pmax(as.numeric(df_model$time_year), tiny_time)
  event <- as.integer(df_model$event_main)
  Z <- model.matrix(stats::as.formula(incidence_formula), data = df_model)
  X <- model.matrix(stats::as.formula(latency_formula), data = df_model)
  
  cure_start <- clamp_prob(km_survival_at(time, event, max(time)))
  gamma_start <- rep(0, ncol(Z))
  names(gamma_start) <- colnames(Z)
  gamma_start[1] <- logit_clip(cure_start)
  
  if (is.null(null_beta_start) || is.null(null_log_aux_start) || anyNA(null_beta_start) || is.na(null_log_aux_start)) {
    null_fit <- fit_no_cure_aft(df_model, latency_formula, family)
    beta_start <- if (isTRUE(null_fit$converged)) null_fit$beta_hat else rep(0, ncol(X))
    log_aux_start <- if (isTRUE(null_fit$converged)) null_fit$log_aux_hat else 0
  } else {
    beta_start <- null_beta_start
    log_aux_start <- null_log_aux_start
  }
  
  gamma_low <- gamma_start
  gamma_low[1] <- logit_clip(0.05)
  gamma_mid <- gamma_start
  gamma_mid[1] <- 0
  
  start_main <- c(gamma_start, beta_start, log_aux_start)
  start_low <- c(gamma_low, beta_start, log_aux_start)
  start_mid <- c(gamma_mid, beta_start, log_aux_start)
  start_zero <- c(rep(0, ncol(Z)), rep(0, ncol(X)), 0)
  
  objective_fn <- function(par) negloglik_mix_aft(par, Z, X, time, event, family)
  
  opt <- run_optim_attempts(
    start_list = list(start_main, start_low, start_mid, start_zero),
    objective_fn = objective_fn,
    control = list(maxit = optim_maxit, reltol = optim_reltol)
  )
  
  if (is.null(opt)) {
    return(list(
      converged = FALSE,
      incidence_formula = incidence_formula,
      latency_formula = latency_formula,
      family = family,
      gamma_hat = NA_real_,
      beta_hat = NA_real_,
      aux_hat = NA_real_,
      log_aux_hat = NA_real_,
      mean_cure_fraction = NA_real_,
      loglik = NA_real_,
      aic = NA_real_,
      optim_par = NA_real_,
      optim_convergence = NA_integer_
    ))
  }
  
  q <- ncol(Z)
  p <- ncol(X)
  gamma_hat <- opt$par[seq_len(q)]
  names(gamma_hat) <- colnames(Z)
  beta_hat <- opt$par[q + seq_len(p)]
  names(beta_hat) <- colnames(X)
  log_aux_hat <- opt$par[q + p + 1L]
  aux_hat <- safe_exp_param(log_aux_hat)[[1]]
  cure_i <- clamp_prob(plogis(drop(Z %*% gamma_hat)))
  loglik <- -opt$value
  
  list(
    converged = isTRUE(opt$convergence == 0L) && is.finite(loglik) && is.finite(aux_hat),
    incidence_formula = incidence_formula,
    latency_formula = latency_formula,
    family = family,
    gamma_hat = gamma_hat,
    beta_hat = beta_hat,
    aux_hat = aux_hat,
    log_aux_hat = log_aux_hat,
    mean_cure_fraction = mean(cure_i),
    loglik = loglik,
    aic = 2 * length(opt$par) - 2 * loglik,
    optim_par = opt$par,
    optim_convergence = opt$convergence
  )
}

get_hsu_formula_spec <- function(dataset_key, formula_registry) {
  if (dataset_key %in% c("PNU", "SNU")) {
    dataset_name <- dataset_key
    latency_rhs <- select_formula_rhs(formula_registry, dataset_name, formula_name = "base")
    incidence_rhs_set <- select_formula_rhs(formula_registry, dataset_name)
    fallback_latency_rhs <- "age_s + sex_num"
    fallback_incidence_rhs <- c("age_s + sex_num", "age_s + sex_num + age_s:sex_num")
    hsu_formula_branch <- "site_free"
    tail_metric_source <- "self"
  } else if (dataset_key == "merged__site_free") {
    dataset_name <- "merged"
    latency_rhs <- select_formula_rhs(formula_registry, dataset_name, formula_name = "base")
    incidence_rhs_set <- select_formula_rhs(formula_registry, dataset_name, site_branch = "site_free")
    fallback_latency_rhs <- "age_s + sex_num"
    fallback_incidence_rhs <- c("age_s + sex_num", "age_s + sex_num + age_s:sex_num")
    hsu_formula_branch <- "site_free"
    tail_metric_source <- "self"
  } else if (dataset_key == "merged__site_adjusted") {
    dataset_name <- "merged"
    latency_rhs <- select_formula_rhs(formula_registry, dataset_name, formula_name = "site_added")
    incidence_rhs_set <- select_formula_rhs(formula_registry, dataset_name, site_branch = "site_adjusted")
    fallback_latency_rhs <- "age_s + sex_num + site"
    fallback_incidence_rhs <- c("age_s + sex_num + site", "age_s + sex_num + age_s:sex_num + site")
    hsu_formula_branch <- "site_adjusted"
    tail_metric_source <- "merged__site_free"
  } else {
    stop(sprintf("Unsupported dataset_key for HSU specification: %s", dataset_key), call. = FALSE)
  }
  
  if (length(latency_rhs) == 0L || all(is.na(latency_rhs))) {
    latency_rhs <- fallback_latency_rhs
  } else {
    latency_rhs <- latency_rhs[[1]]
  }
  
  if (length(incidence_rhs_set) == 0L || all(is.na(incidence_rhs_set))) {
    incidence_rhs_set <- fallback_incidence_rhs
  }
  
  incidence_rhs_set <- limit_formula_set(unique(c("1", incidence_rhs_set)), dataset_key = dataset_key)
  
  list(
    dataset_key = dataset_key,
    source_dataset = dataset_name,
    tail_metric_source = tail_metric_source,
    hsu_formula_branch = hsu_formula_branch,
    latency_formula = paste("~", latency_rhs),
    latency_formula_rhs = latency_rhs,
    incidence_formulas = paste("~", incidence_rhs_set),
    incidence_rhs_set = incidence_rhs_set
  )
}

simulate_from_no_cure_aft <- function(df, null_fit) {
  if (!isTRUE(null_fit$converged)) {
    stop("Cannot simulate bootstrap data from a failed no-cure AFT fit.", call. = FALSE)
  }
  
  df_model <- df
  if ("site" %in% names(df_model)) {
    df_model$site <- factor(df_model$site)
  }
  
  X <- model.matrix(stats::as.formula(null_fit$latency_formula), data = df_model)
  mu <- drop(X %*% null_fit$beta_hat)
  sigma <- null_fit$aux_hat
  family <- null_fit$family
  n <- nrow(df_model)
  
  if (family == "weibull") {
    u <- stats::runif(n, min = 1e-12, max = 1 - 1e-12)
    event_time <- exp(mu) * (-log(u))^sigma
  } else if (family == "lognormal") {
    event_time <- exp(mu + sigma * stats::rnorm(n))
  } else if (family == "loglogistic") {
    u <- stats::runif(n, min = 1e-12, max = 1 - 1e-12)
    event_time <- exp(mu + sigma * stats::qlogis(u))
  } else {
    stop(sprintf("Unsupported simulation family: %s", family), call. = FALSE)
  }
  
  censor_pool <- as.numeric(df$time_year[df$event_main == 0L])
  if (length(censor_pool) >= 2L) {
    censor_time <- sample(censor_pool, size = n, replace = TRUE)
  } else {
    censor_time <- rep(max(as.numeric(df$time_year)) + 0.25, n)
  }
  
  observed_time <- pmin(event_time, censor_time)
  observed_event <- as.integer(event_time <= censor_time)
  
  df_sim <- df
  df_sim$time_year <- observed_time
  df_sim$days_followup <- observed_time * 365.25
  df_sim$event_main <- observed_event
  df_sim$censor_main <- 1L - observed_event
  df_sim$right_censor_flag <- 1L - observed_event
  df_sim$remission_flag <- 0L
  df_sim$status_num <- ifelse(observed_event == 1L, 1L, 0L)
  df_sim$status_label <- factor(ifelse(observed_event == 1L, "transition", "right_censoring"), levels = c("right_censoring", "remission", "transition"))
  
  df_sim
}

compute_hsu_supscore_from_data <- function(df, dataset_key, formula_registry, working_family) {
  spec <- get_hsu_formula_spec(dataset_key, formula_registry)
  null_fit <- fit_no_cure_aft(df, spec$latency_formula, family = working_family)
  
  if (!isTRUE(null_fit$converged)) {
    return(list(
      statistic = NA_real_,
      working_family = working_family,
      selected_formula = NA_character_,
      latency_formula = spec$latency_formula,
      latency_formula_rhs = spec$latency_formula_rhs,
      alt_table = tibble::tibble(),
      null_fit = null_fit,
      alt_fits = list()
    ))
  }
  
  alt_fits <- lapply(spec$incidence_formulas, function(inc_formula) {
    fit_mixture_cure_aft(
      df = df,
      incidence_formula = inc_formula,
      latency_formula = spec$latency_formula,
      family = working_family,
      null_beta_start = null_fit$beta_hat,
      null_log_aux_start = null_fit$log_aux_hat
    )
  })
  
  alt_table <- bind_rows(lapply(seq_along(alt_fits), function(i) {
    fit_obj <- alt_fits[[i]]
    lr_stat <- if (isTRUE(fit_obj$converged)) max(0, 2 * (fit_obj$loglik - null_fit$loglik)) else NA_real_
    
    tibble::tibble(
      working_family = working_family,
      incidence_formula = spec$incidence_formulas[[i]],
      incidence_formula_rhs = strip_formula_rhs(spec$incidence_formulas[[i]]),
      converged = isTRUE(fit_obj$converged),
      loglik = extract_list_scalar(fit_obj, "loglik", NA_real_),
      lr_stat = lr_stat,
      mean_cure_fraction = extract_list_scalar(fit_obj, "mean_cure_fraction", NA_real_)
    )
  }))
  
  lr_for_argmax <- alt_table$lr_stat
  lr_for_argmax[!is.finite(lr_for_argmax)] <- -Inf
  
  if (length(lr_for_argmax) == 0L || all(lr_for_argmax == -Inf)) {
    selected_formula <- NA_character_
    statistic <- 0
  } else {
    selected_index <- which.max(lr_for_argmax)
    selected_formula <- alt_table$incidence_formula[[selected_index]]
    statistic <- max(alt_table$lr_stat, na.rm = TRUE)
  }
  
  list(
    statistic = statistic,
    working_family = working_family,
    selected_formula = selected_formula,
    latency_formula = spec$latency_formula,
    latency_formula_rhs = spec$latency_formula_rhs,
    alt_table = alt_table,
    null_fit = null_fit,
    alt_fits = alt_fits
  )
}

run_hsu_supscore_approx_family <- function(df, dataset_key, formula_registry, bootstrap_reps, seed, working_family) {
  spec <- get_hsu_formula_spec(dataset_key, formula_registry)
  
  if (sum(df$event_main == 1L) < hsu_min_events) {
    return(list(
      statistic = NA_real_,
      p_value = NA_real_,
      adjusted_p_value = NA_real_,
      critical_value = NA_real_,
      selected_formula = NA_character_,
      working_family = working_family,
      latency_formula = spec$latency_formula,
      latency_formula_rhs = spec$latency_formula_rhs,
      alt_table = tibble::tibble(),
      null_fit = NULL,
      bootstrap_stats = rep(NA_real_, bootstrap_reps),
      note = sprintf("HSU approximation skipped because n_event < %s.", hsu_min_events)
    ))
  }
  
  observed <- compute_hsu_supscore_from_data(df, dataset_key, formula_registry, working_family = working_family)
  
  if (!isTRUE(observed$null_fit$converged) || is.na(observed$statistic)) {
    return(list(
      statistic = observed$statistic,
      p_value = NA_real_,
      adjusted_p_value = NA_real_,
      critical_value = NA_real_,
      selected_formula = observed$selected_formula,
      working_family = working_family,
      latency_formula = observed$latency_formula,
      latency_formula_rhs = observed$latency_formula_rhs,
      alt_table = observed$alt_table,
      null_fit = observed$null_fit,
      bootstrap_stats = rep(NA_real_, bootstrap_reps),
      note = sprintf(
        "HSU approximation failed because the no-cure working %s AFT model did not converge.",
        working_family
      )
    ))
  }
  
  label <- paste0("HSU[", dataset_key, ":", working_family, "]")
  bootstrap_stats <- run_bootstrap_chunks(
    n_reps = bootstrap_reps,
    seed = seed,
    label = label,
    worker_fun = function(rep_id) {
      df_boot <- simulate_from_no_cure_aft(df, observed$null_fit)
      compute_hsu_supscore_from_data(df_boot, dataset_key, formula_registry, working_family = working_family)$statistic
    }
  )
  
  n_valid <- sum(is.finite(bootstrap_stats))
  p_value <- if (n_valid > 0L) {
    (sum(bootstrap_stats >= observed$statistic, na.rm = TRUE) + 1) / (n_valid + 1)
  } else {
    NA_real_
  }
  
  critical_value <- if (n_valid > 0L) {
    as.numeric(stats::quantile(bootstrap_stats, probs = 1 - alpha_screening, na.rm = TRUE, type = 8, names = FALSE))
  } else {
    NA_real_
  }
  
  list(
    statistic = observed$statistic,
    p_value = p_value,
    adjusted_p_value = NA_real_,
    critical_value = critical_value,
    selected_formula = observed$selected_formula,
    working_family = working_family,
    latency_formula = observed$latency_formula,
    latency_formula_rhs = observed$latency_formula_rhs,
    alt_table = observed$alt_table,
    null_fit = observed$null_fit,
    bootstrap_stats = bootstrap_stats,
    note = paste0(
      "Self-contained pragmatic approximation to the Hsu sup-score test using a ",
      working_family,
      " AFT working latency model and parametric bootstrap under the no-cure null."
    )
  )
}

run_hsu_familyset_approx <- function(df, dataset_key, formula_registry, bootstrap_reps, seed) {
  family_results <- stats::setNames(
    lapply(seq_along(candidate_hsu_families), function(i) {
      run_hsu_supscore_approx_family(
        df = df,
        dataset_key = dataset_key,
        formula_registry = formula_registry,
        bootstrap_reps = bootstrap_reps,
        seed = seed + i * 1000L,
        working_family = candidate_hsu_families[[i]]
      )
    }),
    candidate_hsu_families
  )
  
  familywise_table <- bind_rows(lapply(candidate_hsu_families, function(family_name) {
    obj <- family_results[[family_name]]
    tibble::tibble(
      working_family = family_name,
      method_name = hsu_method_name(family_name),
      statistic = extract_list_scalar(obj, "statistic", NA_real_),
      raw_p_value = extract_list_scalar(obj, "p_value", NA_real_),
      adjusted_p_value = NA_real_,
      critical_value = extract_list_scalar(obj, "critical_value", NA_real_),
      selected_formula = extract_list_scalar(obj, "selected_formula", NA_character_),
      latency_formula = extract_list_scalar(obj, "latency_formula", NA_character_),
      latency_formula_rhs = extract_list_scalar(obj, "latency_formula_rhs", NA_character_),
      converged = isTRUE(extract_list_scalar(obj$null_fit, "converged", FALSE)),
      note = extract_list_scalar(obj, "note", NA_character_)
    )
  }))
  
  finite_idx <- which(is.finite(familywise_table$raw_p_value))
  if (length(finite_idx) > 0L) {
    familywise_table$adjusted_p_value[finite_idx] <- stats::p.adjust(familywise_table$raw_p_value[finite_idx], method = "holm")
  }
  
  selected_family <- if (length(finite_idx) > 0L) {
    familywise_table$working_family[finite_idx][which.min(familywise_table$raw_p_value[finite_idx])]
  } else {
    NA_character_
  }
  
  selected_formula <- if (!is.na(selected_family)) {
    family_results[[selected_family]]$selected_formula
  } else {
    NA_character_
  }
  
  familyset_p_value <- if (length(finite_idx) > 0L) {
    min(familywise_table$adjusted_p_value[finite_idx], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  finite_stat_idx <- which(is.finite(familywise_table$statistic))
  familyset_statistic <- if (length(finite_stat_idx) > 0L) {
    max(familywise_table$statistic[finite_stat_idx], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  familyset_result <- list(
    statistic = familyset_statistic,
    p_value = familyset_p_value,
    adjusted_p_value = familyset_p_value,
    critical_value = NA_real_,
    selected_family = selected_family,
    selected_formula = selected_formula,
    family_set_name = hsu_family_set_name,
    note = paste0(
      "HSU family-set aggregate based on family-specific pragmatic approximations across ",
      hsu_family_set_name,
      "; official family-set flag uses Holm-adjusted familywise p-values."
    )
  )
  
  list(
    family_results = family_results,
    familywise_table = familywise_table,
    familyset = familyset_result
  )
}

# 🔴 Build: row constructors and decision rules ===============================
classify_presence_p <- function(p_value, alpha = alpha_screening) {
  if (is.na(p_value)) "equivocal" else if (p_value < alpha) "supportive" else "unsupportive"
}

classify_alpha_n <- function(alpha_n, alpha = alpha_screening) {
  if (is.na(alpha_n)) "equivocal" else if (alpha_n < alpha) "supportive" else "unsupportive"
}

classify_xie_flag <- function(p_value, alpha = alpha_screening) {
  if (is.na(p_value)) "equivocal" else if (p_value < alpha) "unsupportive" else "supportive"
}

classify_receus_flag <- function(cure_fraction_hat, receus_ratio_hat) {
  if (is.na(cure_fraction_hat) || is.na(receus_ratio_hat)) {
    return("equivocal")
  }
  if (cure_fraction_hat > receus_pi_threshold && receus_ratio_hat < receus_r_threshold) {
    return("supportive")
  }
  if ((cure_fraction_hat > receus_pi_threshold && receus_ratio_hat < receus_r_equivocal) ||
      (cure_fraction_hat > receus_pi_equivocal && receus_ratio_hat < receus_r_equivocal)) {
    return("equivocal")
  }
  "unsupportive"
}

build_method_row <- function(
    dataset_key,
    source_dataset,
    analysis_variant,
    method_name,
    method_role,
    statistic_name = NA_character_,
    statistic_value = NA_real_,
    p_value = NA_real_,
    adjusted_p_value = NA_real_,
    p_value_type = NA_character_,
    method_flag = "equivocal",
    alpha = alpha_screening,
    model_family = NA_character_,
    selected_family = NA_character_,
    family_set_name = NA_character_,
    latency_formula_rhs = NA_character_,
    incidence_formula_rhs = NA_character_,
    bootstrap_reps = NA_integer_,
    cure_fraction_hat = NA_real_,
    noncure_hat_naive = NA_real_,
    noncure_hat_corrected = NA_real_,
    receus_ratio_hat = NA_real_,
    overall_survival_tau = NA_real_,
    susceptible_survival_tau = NA_real_,
    tau_year = NA_real_,
    aic = NA_real_,
    mz_alpha_n = NA_real_,
    mz_q_n = NA_real_,
    xie_epsilon = NA_real_,
    convergence_flag = NA,
    implementation_label = NA_character_,
    note = NA_character_
) {
  tibble::tibble(
    dataset_key = as.character(dataset_key),
    source_dataset = as.character(source_dataset),
    analysis_variant = as.character(analysis_variant),
    method_name = as.character(method_name),
    method_role = as.character(method_role),
    statistic_name = as.character(statistic_name),
    statistic_value = as.numeric(statistic_value),
    p_value = as.numeric(p_value),
    adjusted_p_value = as.numeric(adjusted_p_value),
    p_value_type = as.character(p_value_type),
    alpha = as.numeric(alpha),
    method_flag = as.character(method_flag),
    model_family = as.character(model_family),
    selected_family = as.character(selected_family),
    family_set_name = as.character(family_set_name),
    latency_formula_rhs = as.character(latency_formula_rhs),
    incidence_formula_rhs = as.character(incidence_formula_rhs),
    bootstrap_reps = as.integer(bootstrap_reps),
    cure_fraction_hat = as.numeric(cure_fraction_hat),
    noncure_hat_naive = as.numeric(noncure_hat_naive),
    noncure_hat_corrected = as.numeric(noncure_hat_corrected),
    receus_ratio_hat = as.numeric(receus_ratio_hat),
    overall_survival_tau = as.numeric(overall_survival_tau),
    susceptible_survival_tau = as.numeric(susceptible_survival_tau),
    tau_year = as.numeric(tau_year),
    aic = as.numeric(aic),
    mz_alpha_n = as.numeric(mz_alpha_n),
    mz_q_n = as.numeric(mz_q_n),
    xie_epsilon = as.numeric(xie_epsilon),
    convergence_flag = as.logical(convergence_flag),
    implementation_label = as.character(implementation_label),
    note = as.character(note)
  )
}

summarise_followup_structure <- function(df, dataset_key, source_dataset, analysis_variant) {
  time <- as.numeric(df$time_year)
  event <- as.integer(df$event_main)
  km_fit <- make_km_fit(time, event)
  
  n_total <- nrow(df)
  n_event <- sum(event == 1L)
  tau <- max(time)
  largest_event_time <- if (n_event > 0L) max(time[event == 1L]) else 0
  
  mz_summary <- compute_mz_extreme_summaries(time, event)
  xie_observed <- compute_xie_statistic(time, event)
  
  tibble::tibble(
    dataset_key = dataset_key,
    source_dataset = source_dataset,
    analysis_variant = analysis_variant,
    screening_context = screening_context,
    n_total = n_total,
    n_event = n_event,
    n_censor_main = sum(df$censor_main == 1L),
    n_right_censor = sum(df$right_censor_flag == 1L),
    n_remission = sum(df$remission_flag == 1L),
    censor_rate_main = safe_divide(sum(df$censor_main == 1L), n_total),
    person_time_year = sum(time),
    tau_year = tau,
    largest_event_time_year = if (n_event > 0L) largest_event_time else NA_real_,
    plateau_length_year = tau - largest_event_time,
    n_censored_after_last_event = if (n_event > 0L) sum(df$censor_main == 1L & time > largest_event_time) else sum(df$censor_main == 1L),
    km_tail_survival = as.numeric(km_survival_at(time, event, tau, km_fit)),
    km_tail_noncure = as.numeric(km_cdf_at(time, event, tau, km_fit)),
    mz_interval_lower_year = mz_summary$lower_bound,
    mz_interval_event_count = mz_summary$N_n,
    mz_alpha_n = mz_summary$alpha_n,
    mz_q_n = mz_summary$q_n,
    xie_epsilon = xie_observed$epsilon,
    xie_p_hat_n = xie_observed$p_hat_n,
    xie_p_hat_G = xie_observed$p_hat_G,
    xie_T_n = xie_observed$T_n,
    xie_degenerate_flag = xie_observed$degenerate,
    summary_note = NA_character_
  )
}

get_method_value <- function(method_df, target_method, target_column) {
  x <- method_df[method_df$method_name == target_method, target_column, drop = TRUE]
  if (length(x) == 0L) NA else x[[1]]
}

derive_dataset_decision_row <- function(dataset_key, source_dataset, analysis_variant, method_df, common_horizons_year, risk_thresholds, screening_note = NA_character_) {
  mz_dn_p <- suppressWarnings(as.numeric(get_method_value(method_df, "maller_zhou_dn", "p_value")))
  hsu_familyset_p <- suppressWarnings(as.numeric(get_method_value(method_df, hsu_familyset_method_name, "p_value")))
  xie_p <- suppressWarnings(as.numeric(get_method_value(method_df, "xie_sufficient_followup", "p_value")))
  receus_weibull_flag <- as.character(get_method_value(method_df, "receus_weibull", "method_flag"))
  receus_aic_flag <- as.character(get_method_value(method_df, "receus_aic", "method_flag"))
  mz_alpha_flag <- as.character(get_method_value(method_df, "maller_zhou_alpha_n", "method_flag"))
  
  presence_p_values <- c(mz_dn_p, hsu_familyset_p)
  presence_p_values <- presence_p_values[!is.na(presence_p_values)]
  presence_support_flag <- length(presence_p_values) > 0L && any(presence_p_values < alpha_screening)
  presence_unsupport_flag <- length(presence_p_values) > 0L && all(presence_p_values >= alpha_screening)
  
  if (is.na(receus_aic_flag) || !nzchar(receus_aic_flag)) {
    receus_aic_flag <- "equivocal"
  }
  
  primary_gate_method <- "RECeUS-AIC"
  primary_gate_flag <- receus_aic_flag
  primary_support_flag <- identical(primary_gate_flag, "supportive")
  primary_unsupport_flag <- identical(primary_gate_flag, "unsupportive")
  primary_equivocal_flag <- identical(primary_gate_flag, "equivocal")
  
  xie_contradiction_flag <- if (is.na(xie_p)) NA else xie_p < alpha_screening
  xie_noncontradiction_flag <- if (is.na(xie_p)) NA else xie_p >= alpha_screening
  
  followup_support_flag <- isTRUE(primary_support_flag) && isTRUE(xie_noncontradiction_flag)
  followup_unsupport_flag <- isTRUE(primary_unsupport_flag) || isTRUE(xie_contradiction_flag)
  
  cure_model_eligibility_flag <- if (primary_support_flag && isTRUE(xie_noncontradiction_flag)) {
    "supportive"
  } else if (primary_unsupport_flag && !presence_support_flag) {
    "unsupportive"
  } else {
    "equivocal"
  }
  
  official_method_names <- c("maller_zhou_dn", "xie_sufficient_followup", "receus_weibull", "receus_aic", hsu_familyset_method_name)
  official_method_df <- method_df[method_df$method_name %in% official_method_names, , drop = FALSE]
  supporting_methods <- compact_character_set(official_method_df$method_name[official_method_df$method_flag == "supportive"])
  contradicting_methods <- compact_character_set(official_method_df$method_name[official_method_df$method_flag == "unsupportive"])
  
  presence_evidence_label <- if (presence_support_flag) {
    "at_least_one_presence_only_test_positive"
  } else if (presence_unsupport_flag) {
    "no_presence_only_test_positive"
  } else {
    "presence_only_evidence_mixed_or_unavailable"
  }
  
  followup_evidence_label <- if (primary_support_flag && isTRUE(xie_noncontradiction_flag)) {
    "receus_aic_supportive_and_not_contradicted_by_xie"
  } else if (primary_support_flag && isTRUE(xie_contradiction_flag)) {
    "receus_aic_supportive_but_contradicted_by_xie"
  } else if (primary_unsupport_flag) {
    "receus_aic_unsupportive"
  } else {
    "receus_aic_equivocal_or_unavailable"
  }
  
  modifier_evidence_label <- if (primary_unsupport_flag && presence_support_flag) {
    "presence_only_signal_upgrades_unsupportive_to_equivocal"
  } else if (primary_support_flag && isTRUE(xie_contradiction_flag)) {
    "xie_contradiction_downgrades_supportive_to_equivocal"
  } else if (primary_equivocal_flag) {
    "primary_gate_equivocal"
  } else {
    "no_modifier_triggered"
  }
  
  primary_gate_flag <- ifelse(is.na(primary_gate_flag) | !nzchar(primary_gate_flag), "equivocal", primary_gate_flag)
  
  tibble::tibble(
    dataset_key = dataset_key,
    source_dataset = source_dataset,
    analysis_variant = analysis_variant,
    screening_context = screening_context,
    decision_rule_version = decision_rule_version,
    primary_gate_method = primary_gate_method,
    primary_gate_flag = primary_gate_flag,
    presence_modifier_flag = as.logical(presence_support_flag),
    followup_contradiction_flag = as.logical(xie_contradiction_flag),
    descriptive_tail_summary_flag = "descriptive_only",
    presence_modifier_label = if (presence_support_flag) "positive_presence_only_signal" else if (presence_unsupport_flag) "no_positive_presence_only_signal" else "mixed_or_unavailable",
    followup_contradiction_label = if (isTRUE(xie_contradiction_flag)) "contradiction_present" else if (isTRUE(xie_noncontradiction_flag)) "no_contradiction" else "unavailable",
    modifier_evidence_label = modifier_evidence_label,
    cure_model_eligibility_flag = cure_model_eligibility_flag,
    final_decision_flag = cure_model_eligibility_flag,
    presence_support_flag = presence_support_flag,
    followup_support_flag = followup_support_flag,
    presence_unsupport_flag = presence_unsupport_flag,
    followup_unsupport_flag = followup_unsupport_flag,
    presence_evidence_label = presence_evidence_label,
    followup_evidence_label = followup_evidence_label,
    maller_zhou_dn_p_value = mz_dn_p,
    hsu_familyset_p_value = hsu_familyset_p,
    xie_centered_bootstrap_p_value = xie_p,
    maller_zhou_alpha_flag = mz_alpha_flag,
    receus_weibull_flag = receus_weibull_flag,
    receus_aic_flag = receus_aic_flag,
    supporting_methods = supporting_methods,
    contradicting_methods = contradicting_methods,
    common_horizon_vector = paste(common_horizons_year, collapse = "|"),
    common_threshold_vector = paste(format(risk_thresholds, trim = TRUE, scientific = FALSE), collapse = "|"),
    carry_forward_stage7 = TRUE,
    carry_forward_stage8 = TRUE,
    screening_repeated_in_stage8 = FALSE,
    screening_note = screening_note
  )
}

build_screening_summary_row <- function(dataset_key, followup_summary, method_results, carry_forward_table) {
  followup_row <- followup_summary[followup_summary$dataset_key == dataset_key, , drop = FALSE]
  method_subset <- method_results[method_results$dataset_key == dataset_key, , drop = FALSE]
  carry_row <- carry_forward_table[carry_forward_table$dataset_key == dataset_key, , drop = FALSE]
  
  tibble::tibble(
    dataset_key = dataset_key,
    source_dataset = extract_frame_cell(carry_row, "source_dataset", NA_character_),
    source_description = extract_frame_cell(carry_row, "source_description", NA_character_),
    analysis_variant = extract_frame_cell(carry_row, "analysis_variant", NA_character_),
    screening_context = extract_frame_cell(carry_row, "screening_context", NA_character_),
    tail_metric_source = extract_frame_cell(carry_row, "tail_metric_source", NA_character_),
    hsu_formula_branch = extract_frame_cell(carry_row, "hsu_formula_branch", NA_character_),
    decision_rule_version = extract_frame_cell(carry_row, "decision_rule_version", NA_character_),
    primary_gate_method = extract_frame_cell(carry_row, "primary_gate_method", NA_character_),
    primary_gate_flag = extract_frame_cell(carry_row, "primary_gate_flag", NA_character_),
    presence_modifier_flag = extract_frame_cell(carry_row, "presence_modifier_flag", NA),
    followup_contradiction_flag = extract_frame_cell(carry_row, "followup_contradiction_flag", NA),
    descriptive_tail_summary_flag = extract_frame_cell(carry_row, "descriptive_tail_summary_flag", NA_character_),
    modifier_evidence_label = extract_frame_cell(carry_row, "modifier_evidence_label", NA_character_),
    n_total = extract_frame_cell(followup_row, "n_total", NA_real_),
    n_event = extract_frame_cell(followup_row, "n_event", NA_real_),
    n_censor_main = extract_frame_cell(followup_row, "n_censor_main", NA_real_),
    n_right_censor = extract_frame_cell(followup_row, "n_right_censor", NA_real_),
    n_remission = extract_frame_cell(followup_row, "n_remission", NA_real_),
    censor_rate_main = extract_frame_cell(followup_row, "censor_rate_main", NA_real_),
    person_time_year = extract_frame_cell(followup_row, "person_time_year", NA_real_),
    tau_year = extract_frame_cell(followup_row, "tau_year", NA_real_),
    largest_event_time_year = extract_frame_cell(followup_row, "largest_event_time_year", NA_real_),
    plateau_length_year = extract_frame_cell(followup_row, "plateau_length_year", NA_real_),
    n_censored_after_last_event = extract_frame_cell(followup_row, "n_censored_after_last_event", NA_real_),
    km_tail_survival = extract_frame_cell(followup_row, "km_tail_survival", NA_real_),
    km_tail_noncure = extract_frame_cell(followup_row, "km_tail_noncure", NA_real_),
    mz_interval_lower_year = extract_frame_cell(followup_row, "mz_interval_lower_year", NA_real_),
    mz_interval_event_count = extract_frame_cell(followup_row, "mz_interval_event_count", NA_real_),
    maller_zhou_alpha_n = get_method_value(method_subset, "maller_zhou_alpha_n", "statistic_value"),
    maller_zhou_alpha_flag = get_method_value(method_subset, "maller_zhou_alpha_n", "method_flag"),
    maller_zhou_q_n = extract_frame_cell(followup_row, "mz_q_n", NA_real_),
    maller_zhou_dn_stat = get_method_value(method_subset, "maller_zhou_dn", "statistic_value"),
    maller_zhou_dn_p_value = get_method_value(method_subset, "maller_zhou_dn", "p_value"),
    maller_zhou_dn_flag = get_method_value(method_subset, "maller_zhou_dn", "method_flag"),
    xie_epsilon = extract_frame_cell(followup_row, "xie_epsilon", NA_real_),
    xie_p_hat_n = extract_frame_cell(followup_row, "xie_p_hat_n", NA_real_),
    xie_p_hat_G = extract_frame_cell(followup_row, "xie_p_hat_G", NA_real_),
    xie_T_n = get_method_value(method_subset, "xie_sufficient_followup", "statistic_value"),
    xie_centered_bootstrap_p_value = get_method_value(method_subset, "xie_sufficient_followup", "p_value"),
    xie_flag = get_method_value(method_subset, "xie_sufficient_followup", "method_flag"),
    receus_weibull_family = get_method_value(method_subset, "receus_weibull", "selected_family"),
    receus_weibull_aic = get_method_value(method_subset, "receus_weibull", "aic"),
    receus_weibull_cure_fraction = get_method_value(method_subset, "receus_weibull", "cure_fraction_hat"),
    receus_weibull_ratio = get_method_value(method_subset, "receus_weibull", "receus_ratio_hat"),
    receus_weibull_flag = get_method_value(method_subset, "receus_weibull", "method_flag"),
    receus_aic_selected_family = get_method_value(method_subset, "receus_aic", "selected_family"),
    receus_aic_aic = get_method_value(method_subset, "receus_aic", "aic"),
    receus_aic_cure_fraction = get_method_value(method_subset, "receus_aic", "cure_fraction_hat"),
    receus_aic_ratio = get_method_value(method_subset, "receus_aic", "receus_ratio_hat"),
    receus_aic_flag = get_method_value(method_subset, "receus_aic", "method_flag"),
    hsu_weibull_stat = get_method_value(method_subset, "hsu_supscore_weibull_approx", "statistic_value"),
    hsu_weibull_p_value = get_method_value(method_subset, "hsu_supscore_weibull_approx", "p_value"),
    hsu_weibull_adjusted_p_value = get_method_value(method_subset, "hsu_supscore_weibull_approx", "adjusted_p_value"),
    hsu_weibull_flag = get_method_value(method_subset, "hsu_supscore_weibull_approx", "method_flag"),
    hsu_lognormal_stat = get_method_value(method_subset, "hsu_supscore_lognormal_approx", "statistic_value"),
    hsu_lognormal_p_value = get_method_value(method_subset, "hsu_supscore_lognormal_approx", "p_value"),
    hsu_lognormal_adjusted_p_value = get_method_value(method_subset, "hsu_supscore_lognormal_approx", "adjusted_p_value"),
    hsu_lognormal_flag = get_method_value(method_subset, "hsu_supscore_lognormal_approx", "method_flag"),
    hsu_loglogistic_stat = get_method_value(method_subset, "hsu_supscore_loglogistic_approx", "statistic_value"),
    hsu_loglogistic_p_value = get_method_value(method_subset, "hsu_supscore_loglogistic_approx", "p_value"),
    hsu_loglogistic_adjusted_p_value = get_method_value(method_subset, "hsu_supscore_loglogistic_approx", "adjusted_p_value"),
    hsu_loglogistic_flag = get_method_value(method_subset, "hsu_supscore_loglogistic_approx", "method_flag"),
    hsu_familyset_stat = get_method_value(method_subset, hsu_familyset_method_name, "statistic_value"),
    hsu_familyset_p_value = get_method_value(method_subset, hsu_familyset_method_name, "p_value"),
    hsu_familyset_flag = get_method_value(method_subset, hsu_familyset_method_name, "method_flag"),
    hsu_familyset_name = get_method_value(method_subset, hsu_familyset_method_name, "family_set_name"),
    hsu_familyset_selected_family = get_method_value(method_subset, hsu_familyset_method_name, "selected_family"),
    hsu_latency_formula_rhs = get_method_value(method_subset, hsu_familyset_method_name, "latency_formula_rhs"),
    hsu_selected_incidence_formula_rhs = get_method_value(method_subset, hsu_familyset_method_name, "incidence_formula_rhs"),
    cure_model_eligibility_flag = extract_frame_cell(carry_row, "cure_model_eligibility_flag", NA_character_),
    final_decision_flag = extract_frame_cell(carry_row, "final_decision_flag", NA_character_),
    presence_evidence_label = extract_frame_cell(carry_row, "presence_evidence_label", NA_character_),
    followup_evidence_label = extract_frame_cell(carry_row, "followup_evidence_label", NA_character_),
    supporting_methods = extract_frame_cell(carry_row, "supporting_methods", NA_character_),
    contradicting_methods = extract_frame_cell(carry_row, "contradicting_methods", NA_character_),
    followup_summary_note = extract_frame_cell(followup_row, "summary_note", NA_character_),
    screening_note = extract_frame_cell(carry_row, "screening_note", NA_character_)
  )
}

# 🔴 Define: sharded-bootstrap orchestration helpers ===============================
build_variant_registry <- function(source_lookup, formula_registry) {
  variant_specs <- stats::setNames(
    lapply(dataset_order_full, function(k) get_hsu_formula_spec(k, formula_registry)),
    dataset_order_full
  )
  
  tibble::tibble(
    dataset_key = dataset_order_full,
    source_dataset = c("PNU", "SNU", "merged", "merged"),
    source_description = c(
      unname(source_lookup[["PNU"]]),
      unname(source_lookup[["SNU"]]),
      unname(source_lookup[["merged"]]),
      unname(source_lookup[["merged"]])
    ),
    analysis_variant = c("single_cohort", "single_cohort", "site_free", "site_adjusted"),
    screening_context = screening_context,
    stage1_dataset_name = c("PNU", "SNU", "merged", "merged"),
    tail_metric_source = vapply(variant_specs, function(x) x$tail_metric_source, character(1)),
    hsu_formula_branch = vapply(variant_specs, function(x) x$hsu_formula_branch, character(1)),
    hsu_latency_formula_rhs = vapply(variant_specs, function(x) x$latency_formula_rhs, character(1)),
    hsu_incidence_formula_rhs_set = vapply(variant_specs, function(x) compact_character_set(x$incidence_rhs_set), character(1)),
    hsu_incidence_formula_count = vapply(variant_specs, function(x) length(x$incidence_rhs_set), integer(1)),
    hsu_working_family_set = hsu_family_set_name,
    screening_variant_note = c(
      "Direct single-cohort screening on Stage 1 PNU data.",
      "Direct single-cohort screening on Stage 1 SNU data.",
      "Direct merged-data screening for the site-free structural branch.",
      "Raw-tail metrics inherit the merged site-free data; only the HSU working-model branch is rerun with site adjustment."
    )
  ) %>%
    arrange(match(dataset_key, dataset_order_full))
}

build_base_variant_map <- function(prepared_datasets) {
  list(
    PNU = list(data = prepared_datasets[["PNU"]], source_dataset = "PNU", analysis_variant = "single_cohort"),
    SNU = list(data = prepared_datasets[["SNU"]], source_dataset = "SNU", analysis_variant = "single_cohort"),
    merged__site_free = list(data = prepared_datasets[["merged"]], source_dataset = "merged", analysis_variant = "site_free")
  )
}

build_dataset_map_all <- function(prepared_datasets) {
  list(
    PNU = list(data = prepared_datasets[["PNU"]], source_dataset = "PNU", analysis_variant = "single_cohort"),
    SNU = list(data = prepared_datasets[["SNU"]], source_dataset = "SNU", analysis_variant = "single_cohort"),
    merged__site_free = list(data = prepared_datasets[["merged"]], source_dataset = "merged", analysis_variant = "site_free"),
    merged__site_adjusted = list(data = prepared_datasets[["merged"]], source_dataset = "merged", analysis_variant = "site_adjusted")
  )
}

build_shard_registry <- function() {
  registry_rows <- list()
  created_now <- now_utc_chr()
  
  xie_keys <- c("PNU", "SNU", "merged__site_free")
  xie_ranges <- split_rep_ranges(bootstrap_reps_xie, xie_shard_size)
  
  if (length(xie_ranges) > 0L) {
    for (dataset_key in xie_keys) {
      for (shard_idx in seq_along(xie_ranges)) {
        rep_ids <- xie_ranges[[shard_idx]]
        registry_rows[[length(registry_rows) + 1L]] <- tibble::tibble(
          task_id = sprintf("xie__%s__%03d", dataset_key, shard_idx),
          task_type = "xie",
          dataset_key = dataset_key,
          source_dataset = source_dataset_from_key(dataset_key),
          analysis_variant = analysis_variant_from_key(dataset_key),
          working_family = NA_character_,
          shard_id = as.character(shard_idx),
          n_shards = as.character(length(xie_ranges)),
          rep_id_start = as.character(min(rep_ids)),
          rep_id_end = as.character(max(rep_ids)),
          n_reps = as.character(length(rep_ids)),
          status = "pending",
          worker_id = NA_character_,
          attempt_count = "0",
          created_at_utc = created_now,
          started_at_utc = NA_character_,
          heartbeat_at_utc = NA_character_,
          finished_at_utc = NA_character_,
          n_valid_stats = NA_character_,
          shard_file = make_xie_shard_file(dataset_key, shard_idx, length(xie_ranges)),
          observed_file = make_xie_observed_file(dataset_key),
          note = "Xie bootstrap shard"
        )
      }
    }
  }
  
  hsu_ranges <- split_rep_ranges(bootstrap_reps_hsu, hsu_shard_size)
  if (length(hsu_ranges) > 0L) {
    for (dataset_key in dataset_order_full) {
      for (family_name in candidate_hsu_families) {
        for (shard_idx in seq_along(hsu_ranges)) {
          rep_ids <- hsu_ranges[[shard_idx]]
          registry_rows[[length(registry_rows) + 1L]] <- tibble::tibble(
            task_id = sprintf("hsu__%s__%s__%03d", dataset_key, family_name, shard_idx),
            task_type = "hsu",
            dataset_key = dataset_key,
            source_dataset = source_dataset_from_key(dataset_key),
            analysis_variant = analysis_variant_from_key(dataset_key),
            working_family = family_name,
            shard_id = as.character(shard_idx),
            n_shards = as.character(length(hsu_ranges)),
            rep_id_start = as.character(min(rep_ids)),
            rep_id_end = as.character(max(rep_ids)),
            n_reps = as.character(length(rep_ids)),
            status = "pending",
            worker_id = NA_character_,
            attempt_count = "0",
            created_at_utc = created_now,
            started_at_utc = NA_character_,
            heartbeat_at_utc = NA_character_,
            finished_at_utc = NA_character_,
            n_valid_stats = NA_character_,
            shard_file = make_hsu_shard_file(dataset_key, family_name, shard_idx, length(hsu_ranges)),
            observed_file = make_hsu_observed_file(dataset_key, family_name),
            note = "HSU bootstrap shard"
          )
        }
      }
    }
  }
  
  normalize_shard_registry(compact_bind_rows(registry_rows)) %>%
    arrange(task_type, match(dataset_key, dataset_order_full), working_family, as_integer_or_na(shard_id))
}

compute_xie_observed_object <- function(df, dataset_key, source_dataset, analysis_variant, seed_base) {
  list(
    dataset_key = dataset_key,
    source_dataset = source_dataset,
    analysis_variant = analysis_variant,
    seed_base = as.integer(seed_base),
    tau_year = max(as.numeric(df$time_year)),
    time = as.numeric(df$time_year),
    event = as.integer(df$event_main),
    observed = compute_xie_statistic(as.numeric(df$time_year), as.integer(df$event_main))
  )
}

finalize_xie_from_bootstrap_stats <- function(observed_object, bootstrap_stats) {
  observed <- observed_object$observed
  bootstrap_stats <- as.numeric(bootstrap_stats)
  
  if (is.na(observed$T_n)) {
    return(list(
      observed = observed,
      p_centered = NA_real_,
      p_raw = NA_real_,
      critical_centered = NA_real_,
      bootstrap_stats = bootstrap_stats
    ))
  }
  
  if (isTRUE(observed$degenerate) && observed$T_n == 0) {
    return(list(
      observed = observed,
      p_centered = 1,
      p_raw = 1,
      critical_centered = 0,
      bootstrap_stats = if (length(bootstrap_stats) == 0L) rep(0, bootstrap_reps_xie) else bootstrap_stats
    ))
  }
  
  centered_stats <- bootstrap_stats - observed$T_n
  n_valid_centered <- sum(is.finite(centered_stats))
  n_valid_raw <- sum(is.finite(bootstrap_stats))
  
  critical_centered <- if (n_valid_centered > 0L) {
    as.numeric(stats::quantile(centered_stats, probs = 1 - alpha_screening, na.rm = TRUE, type = 8, names = FALSE))
  } else {
    NA_real_
  }
  
  p_centered <- if (n_valid_centered > 0L) {
    (sum(centered_stats >= observed$T_n, na.rm = TRUE) + 1) / (n_valid_centered + 1)
  } else {
    NA_real_
  }
  
  p_raw <- if (n_valid_raw > 0L) {
    (sum(bootstrap_stats >= observed$T_n, na.rm = TRUE) + 1) / (n_valid_raw + 1)
  } else {
    NA_real_
  }
  
  list(
    observed = observed,
    p_centered = p_centered,
    p_raw = p_raw,
    critical_centered = critical_centered,
    bootstrap_stats = bootstrap_stats
  )
}

finalize_hsu_family_from_bootstrap_stats <- function(observed_object, bootstrap_stats) {
  observed <- observed_object$result
  bootstrap_stats <- as.numeric(bootstrap_stats)
  working_family <- observed_object$working_family
  
  if (sum(is.finite(bootstrap_stats)) == 0L || is.null(observed$null_fit) || !isTRUE(observed$null_fit$converged) || is.na(observed$statistic)) {
    return(list(
      statistic = observed$statistic,
      p_value = NA_real_,
      adjusted_p_value = NA_real_,
      critical_value = NA_real_,
      selected_formula = observed$selected_formula,
      working_family = working_family,
      latency_formula = observed$latency_formula,
      latency_formula_rhs = observed$latency_formula_rhs,
      alt_table = observed$alt_table,
      null_fit = observed$null_fit,
      bootstrap_stats = bootstrap_stats,
      note = sprintf(
        "HSU approximation failed or bootstrap statistics were unavailable for the %s AFT working family.",
        working_family
      )
    ))
  }
  
  n_valid <- sum(is.finite(bootstrap_stats))
  p_value <- if (n_valid > 0L) {
    (sum(bootstrap_stats >= observed$statistic, na.rm = TRUE) + 1) / (n_valid + 1)
  } else {
    NA_real_
  }
  
  critical_value <- if (n_valid > 0L) {
    as.numeric(stats::quantile(bootstrap_stats, probs = 1 - alpha_screening, na.rm = TRUE, type = 8, names = FALSE))
  } else {
    NA_real_
  }
  
  list(
    statistic = observed$statistic,
    p_value = p_value,
    adjusted_p_value = NA_real_,
    critical_value = critical_value,
    selected_formula = observed$selected_formula,
    working_family = working_family,
    latency_formula = observed$latency_formula,
    latency_formula_rhs = observed$latency_formula_rhs,
    alt_table = observed$alt_table,
    null_fit = observed$null_fit,
    bootstrap_stats = bootstrap_stats,
    note = paste0(
      "Self-contained pragmatic approximation to the Hsu sup-score test using a ",
      working_family,
      " AFT working latency model and sharded parametric bootstrap under the no-cure null."
    )
  )
}

reduce_hsu_familyset_from_family_results <- function(family_results) {
  familywise_table <- bind_rows(lapply(candidate_hsu_families, function(family_name) {
    obj <- family_results[[family_name]]
    tibble::tibble(
      working_family = family_name,
      method_name = hsu_method_name(family_name),
      statistic = extract_list_scalar(obj, "statistic", NA_real_),
      raw_p_value = extract_list_scalar(obj, "p_value", NA_real_),
      adjusted_p_value = NA_real_,
      critical_value = extract_list_scalar(obj, "critical_value", NA_real_),
      selected_formula = extract_list_scalar(obj, "selected_formula", NA_character_),
      latency_formula = extract_list_scalar(obj, "latency_formula", NA_character_),
      latency_formula_rhs = extract_list_scalar(obj, "latency_formula_rhs", NA_character_),
      converged = isTRUE(extract_list_scalar(obj$null_fit, "converged", FALSE)),
      note = extract_list_scalar(obj, "note", NA_character_)
    )
  }))
  
  finite_idx <- which(is.finite(familywise_table$raw_p_value))
  if (length(finite_idx) > 0L) {
    familywise_table$adjusted_p_value[finite_idx] <- stats::p.adjust(familywise_table$raw_p_value[finite_idx], method = "holm")
  }
  
  selected_family <- if (length(finite_idx) > 0L) {
    familywise_table$working_family[finite_idx][which.min(familywise_table$raw_p_value[finite_idx])]
  } else {
    NA_character_
  }
  
  selected_formula <- if (!is.na(selected_family)) {
    family_results[[selected_family]]$selected_formula
  } else {
    NA_character_
  }
  
  familyset_p_value <- if (length(finite_idx) > 0L) {
    min(familywise_table$adjusted_p_value[finite_idx], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  finite_stat_idx <- which(is.finite(familywise_table$statistic))
  familyset_statistic <- if (length(finite_stat_idx) > 0L) {
    max(familywise_table$statistic[finite_stat_idx], na.rm = TRUE)
  } else {
    NA_real_
  }
  
  familyset_result <- list(
    statistic = familyset_statistic,
    p_value = familyset_p_value,
    adjusted_p_value = familyset_p_value,
    critical_value = NA_real_,
    selected_family = selected_family,
    selected_formula = selected_formula,
    family_set_name = hsu_family_set_name,
    note = paste0(
      "HSU family-set aggregate based on family-specific sharded pragmatic approximations across ",
      hsu_family_set_name,
      "; official family-set flag uses Holm-adjusted familywise p-values."
    )
  )
  
  list(
    familywise_table = familywise_table,
    familyset = familyset_result
  )
}

collect_bootstrap_stats_from_registry <- function(registry_df, task_type, dataset_key, total_reps, working_family = NA_character_) {
  subset_df <- registry_df %>%
    filter(
      task_type == !!task_type,
      dataset_key == !!dataset_key,
      status == "completed"
    )
  
  if (!is.na(working_family)) {
    subset_df <- subset_df %>% filter(working_family == !!working_family)
  }
  
  out <- rep(NA_real_, total_reps)
  if (nrow(subset_df) == 0L) {
    return(out)
  }
  
  for (ii in seq_len(nrow(subset_df))) {
    shard_file <- subset_df$shard_file[[ii]]
    if (!file.exists(shard_file)) {
      next
    }
    obj <- readRDS(shard_file)
    rep_ids <- as_integer_or_na(obj$rep_ids)
    stats <- as.numeric(obj$statistics)
    valid <- which(!is.na(rep_ids) & rep_ids >= 1L & rep_ids <= total_reps)
    if (length(valid) > 0L) {
      out[rep_ids[valid]] <- stats[valid]
    }
  }
  
  out
}

requeue_stale_and_retryable_tasks <- function(registry_df) {
  registry_df <- normalize_shard_registry(registry_df)
  
  running_idx <- which(registry_df$status == "running")
  if (length(running_idx) > 0L) {
    for (idx in running_idx) {
      worker_for_row <- registry_df$worker_id[[idx]]
      fallback_ts <- coalesce_chr_scalar(registry_df$heartbeat_at_utc[[idx]], registry_df$started_at_utc[[idx]])
      last_seen_age <- worker_last_seen_age_sec(worker_for_row, fallback_timestamp = fallback_ts)
      if (is.finite(last_seen_age) && last_seen_age > task_stale_after_sec) {
        registry_df$status[[idx]] <- "pending"
        registry_df$worker_id[[idx]] <- NA_character_
        registry_df$started_at_utc[[idx]] <- NA_character_
        registry_df$heartbeat_at_utc[[idx]] <- NA_character_
        registry_df$finished_at_utc[[idx]] <- NA_character_
        registry_df$n_valid_stats[[idx]] <- NA_character_
        registry_df$note[[idx]] <- append_note(
          registry_df$note[[idx]],
          paste0(
            "Automatically requeued stale running shard from worker `",
            worker_for_row,
            "` after ",
            round(last_seen_age),
            " seconds without heartbeat."
          )
        )
      }
    }
  }
  
  failed_idx <- which(registry_df$status == "failed")
  if (length(failed_idx) > 0L) {
    for (idx in failed_idx) {
      attempts <- as_integer_or_na(registry_df$attempt_count[[idx]])
      attempts <- ifelse(is.na(attempts), 0L, attempts)
      fail_age <- timestamp_age_sec(registry_df$finished_at_utc[[idx]], default = Inf)
      if (attempts < max_task_attempts && is.finite(fail_age) && fail_age >= failed_task_retry_after_sec) {
        registry_df$status[[idx]] <- "pending"
        registry_df$worker_id[[idx]] <- NA_character_
        registry_df$started_at_utc[[idx]] <- NA_character_
        registry_df$heartbeat_at_utc[[idx]] <- NA_character_
        registry_df$finished_at_utc[[idx]] <- NA_character_
        registry_df$n_valid_stats[[idx]] <- NA_character_
        registry_df$note[[idx]] <- append_note(
          registry_df$note[[idx]],
          paste0(
            "Automatically requeued failed shard for retry ",
            attempts + 1L,
            " of ",
            max_task_attempts,
            "."
          )
        )
      }
    }
  }
  
  registry_df
}

claim_next_task <- function() {
  with_dir_lock(registry_lock_dir, {
    if (!file.exists(shard_registry_file)) {
      return(NULL)
    }
    registry_df <- requeue_stale_and_retryable_tasks(read_shard_registry())
    pending_idx <- which(registry_df$status == "pending")
    if (length(pending_idx) == 0L) {
      write_shard_registry(registry_df)
      return(NULL)
    }
    idx <- pending_idx[[1]]
    now_chr <- now_utc_chr()
    attempts <- as_integer_or_na(registry_df$attempt_count[[idx]])
    attempts <- ifelse(is.na(attempts), 0L, attempts) + 1L
    registry_df$status[[idx]] <- "running"
    registry_df$worker_id[[idx]] <- worker_id
    registry_df$attempt_count[[idx]] <- as.character(attempts)
    registry_df$started_at_utc[[idx]] <- now_chr
    registry_df$heartbeat_at_utc[[idx]] <- now_chr
    registry_df$finished_at_utc[[idx]] <- NA_character_
    write_shard_registry(registry_df)
    task_row <- registry_df[idx, , drop = FALSE]
    touch_worker_heartbeat(
      task_id = task_row$task_id[[1]],
      dataset_key = task_row$dataset_key[[1]],
      stage = "claimed",
      note = paste("Claimed", task_row$task_id[[1]])
    )
    task_row
  }, timeout_sec = worker_claim_lock_timeout_sec, poll_sec = 0.25)
}

update_task_status <- function(task_id, status, n_valid_stats = NA_integer_, note = NA_character_) {
  with_dir_lock(registry_lock_dir, {
    registry_df <- read_shard_registry()
    idx <- which(registry_df$task_id == task_id)
    if (length(idx) == 0L) {
      stop(sprintf("Task id not found in shard registry: %s", task_id), call. = FALSE)
    }
    idx <- idx[[1]]
    now_chr <- now_utc_chr()
    registry_df$status[[idx]] <- status
    registry_df$worker_id[[idx]] <- worker_id
    registry_df$heartbeat_at_utc[[idx]] <- now_chr
    if (!is.na(n_valid_stats)) {
      registry_df$n_valid_stats[[idx]] <- as.character(as.integer(n_valid_stats))
    }
    if (!is.na(note) && nzchar(note)) {
      registry_df$note[[idx]] <- append_note(registry_df$note[[idx]], note)
    }
    if (identical(status, "running")) {
      registry_df$started_at_utc[[idx]] <- now_chr
    } else {
      registry_df$finished_at_utc[[idx]] <- now_chr
    }
    write_shard_registry(registry_df)
  }, timeout_sec = 600, poll_sec = 0.25)
  invisible(TRUE)
}

run_xie_shard_task <- function(task_row) {
  observed_object <- readRDS(task_row$observed_file[[1]])
  rep_ids <- seq.int(as_integer_or_na(task_row$rep_id_start[[1]]), as_integer_or_na(task_row$rep_id_end[[1]]))
  n <- length(observed_object$time)
  statistics <- rep(NA_real_, length(rep_ids))
  
  for (ii in seq_along(rep_ids)) {
    rep_id <- rep_ids[[ii]]
    if (ii == 1L || ii %% task_heartbeat_every_reps == 0L) {
      touch_worker_heartbeat(
        task_id = task_row$task_id[[1]],
        dataset_key = task_row$dataset_key[[1]],
        stage = "running",
        note = sprintf("Xie shard progress %s rep %s/%s", task_row$task_id[[1]], ii, length(rep_ids))
      )
    }
    set.seed(observed_object$seed_base + rep_id)
    draw_index <- sample.int(n = n, size = n, replace = TRUE)
    statistics[[ii]] <- compute_xie_statistic(
      observed_object$time[draw_index],
      observed_object$event[draw_index]
    )$T_n
  }
  
  list(
    task_id = task_row$task_id[[1]],
    task_type = "xie",
    dataset_key = task_row$dataset_key[[1]],
    working_family = NA_character_,
    rep_ids = rep_ids,
    statistics = statistics
  )
}

run_hsu_shard_task <- function(task_row, state_bundle) {
  observed_object <- readRDS(task_row$observed_file[[1]])
  prepared_datasets <- readRDS(prepared_datasets_file)
  source_dataset <- observed_object$source_dataset
  df <- prepared_datasets[[source_dataset]]
  rep_ids <- seq.int(as_integer_or_na(task_row$rep_id_start[[1]]), as_integer_or_na(task_row$rep_id_end[[1]]))
  formula_registry <- state_bundle$inputs$formula_registry
  statistics <- rep(NA_real_, length(rep_ids))
  
  for (ii in seq_along(rep_ids)) {
    rep_id <- rep_ids[[ii]]
    if (ii == 1L || ii %% task_heartbeat_every_reps == 0L) {
      touch_worker_heartbeat(
        task_id = task_row$task_id[[1]],
        dataset_key = task_row$dataset_key[[1]],
        stage = "running",
        note = sprintf("HSU shard progress %s rep %s/%s", task_row$task_id[[1]], ii, length(rep_ids))
      )
    }
    set.seed(observed_object$seed_base + rep_id)
    df_boot <- simulate_from_no_cure_aft(df, observed_object$result$null_fit)
    statistics[[ii]] <- compute_hsu_supscore_from_data(
      df = df_boot,
      dataset_key = observed_object$dataset_key,
      formula_registry = formula_registry,
      working_family = observed_object$working_family
    )$statistic
  }
  
  list(
    task_id = task_row$task_id[[1]],
    task_type = "hsu",
    dataset_key = task_row$dataset_key[[1]],
    working_family = task_row$working_family[[1]],
    rep_ids = rep_ids,
    statistics = statistics
  )
}

run_one_claimed_task <- function(task_row) {
  task_type <- task_row$task_type[[1]]
  task_id <- task_row$task_id[[1]]
  state_bundle <- readRDS(state_bundle_file)
  log_step("Claimed shard: ", task_id, " (", task_type, ")")
  write_runtime_status("task_started", dataset_key = task_row$dataset_key[[1]], note = paste("Claimed", task_id))
  
  result <- tryCatch(
    {
      if (identical(task_type, "xie")) {
        run_xie_shard_task(task_row)
      } else if (identical(task_type, "hsu")) {
        run_hsu_shard_task(task_row, state_bundle = state_bundle)
      } else {
        stop(sprintf("Unsupported task_type: %s", task_type), call. = FALSE)
      }
    },
    error = function(e) {
      update_task_status(task_id, status = "failed", n_valid_stats = 0L, note = conditionMessage(e))
      touch_worker_heartbeat(
        task_id = task_id,
        dataset_key = task_row$dataset_key[[1]],
        stage = "failed",
        note = conditionMessage(e)
      )
      log_step("Task failed: ", task_id, " :: ", conditionMessage(e), level = "WARN")
      write_runtime_status("task_failed", dataset_key = task_row$dataset_key[[1]], note = paste(task_id, conditionMessage(e)))
      return(NULL)
    }
  )
  
  if (is.null(result)) {
    return(invisible(FALSE))
  }
  
  safe_save_rds_atomic(result, task_row$shard_file[[1]])
  update_task_status(task_id, status = "completed", n_valid_stats = sum(is.finite(result$statistics)))
  touch_worker_heartbeat(
    task_id = task_id,
    dataset_key = task_row$dataset_key[[1]],
    stage = "completed",
    note = paste("Completed", task_id)
  )
  log_step(
    "Completed shard: ", task_id,
    " with ", sum(is.finite(result$statistics)), "/", length(result$statistics), " valid statistics."
  )
  write_runtime_status("task_completed", dataset_key = task_row$dataset_key[[1]], note = paste("Completed", task_id))
  invisible(TRUE)
}

prepare_observed_and_registry <- function() {
  ensure_export_path_writable(export_path)
  
  if (isTRUE(refresh_run_state)) {
    reset_run_state_files()
  }
  
  if (file.exists(runtime_status_file)) {
    unlink(runtime_status_file)
  }
  writeLines(character(0), progress_log_file)
  dir.create(heartbeat_dir, recursive = TRUE, showWarnings = FALSE)
  existing_heartbeat_files <- list.files(heartbeat_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (length(existing_heartbeat_files) > 0L) {
    unlink(existing_heartbeat_files, recursive = TRUE, force = TRUE)
  }
  
  log_step("Stage 6 sharded run initialized. run_mode=", run_mode)
  log_step(
    "Controls: xie_bootstrap=", bootstrap_reps_xie,
    ", hsu_bootstrap=", bootstrap_reps_hsu,
    ", xie_shard_size=", xie_shard_size,
    ", hsu_shard_size=", hsu_shard_size,
    ", parallel_cores=", parallel_cores,
    ", hsu_families=", hsu_family_set_name,
    ", task_stale_after_sec=", task_stale_after_sec,
    ", max_task_attempts=", max_task_attempts
  )
  write_runtime_status("initialized", note = "Stage 6 sharded run initialized.")
  
  ingest <- ingest_stage1_assets()
  
  prepared_datasets <- list(
    PNU = prepare_stage6_dataset(ingest$analysis_datasets_stage1[["PNU"]], "PNU", ingest$scaling_registry),
    SNU = prepare_stage6_dataset(ingest$analysis_datasets_stage1[["SNU"]], "SNU", ingest$scaling_registry),
    merged = prepare_stage6_dataset(ingest$analysis_datasets_stage1[["merged"]], "merged", ingest$scaling_registry)
  )
  
  safe_save_rds_atomic(prepared_datasets, prepared_datasets_file)
  
  variant_registry <- build_variant_registry(
    source_lookup = ingest$source_lookup,
    formula_registry = ingest$formula_registry
  )
  
  base_variant_map <- build_base_variant_map(prepared_datasets)
  all_variant_map <- build_dataset_map_all(prepared_datasets)
  
  followup_summary_list <- list()
  nonbootstrap_method_row_list <- list()
  receus_candidate_list <- list()
  
  fit_objects_base <- list(
    maller_zhou = list(),
    receus = list(),
    xie_observed = list(),
    hsu_observed = list(),
    inheritance_notes = list()
  )
  
  for (i in seq_along(base_variant_map)) {
    current_key <- names(base_variant_map)[[i]]
    current_map <- base_variant_map[[i]]
    current_data <- current_map$data
    current_source <- current_map$source_dataset
    current_variant <- current_map$analysis_variant
    current_tau <- max(as.numeric(current_data$time_year))
    
    log_step("Observed-only setup: ", current_key)
    
    followup_summary_list[[length(followup_summary_list) + 1L]] <- summarise_followup_structure(
      df = current_data,
      dataset_key = current_key,
      source_dataset = current_source,
      analysis_variant = current_variant
    )
    
    mz_object <- run_maller_zhou_presence(current_data)
    fit_objects_base$maller_zhou[[current_key]] <- mz_object
    
    nonbootstrap_method_row_list[[length(nonbootstrap_method_row_list) + 1L]] <- build_method_row(
      dataset_key = current_key,
      source_dataset = current_source,
      analysis_variant = current_variant,
      method_name = "maller_zhou_dn",
      method_role = "presence_only",
      statistic_name = "d_n",
      statistic_value = mz_object$result$dn,
      p_value = mz_object$result$p_value,
      p_value_type = "boundary_mixture_chisq",
      method_flag = classify_presence_p(mz_object$result$p_value),
      alpha = alpha_screening,
      model_family = "exponential",
      selected_family = "exponential",
      tau_year = current_tau,
      cure_fraction_hat = mz_object$result$cure_fraction_hat,
      convergence_flag = isTRUE(extract_list_scalar(mz_object$fit_object$alt, "converged", FALSE)),
      implementation_label = "exact_parametric_boundary_test",
      note = "Parametric Maller-Zhou boundary deviance test under an exponential mixture-cure working model."
    )
    
    nonbootstrap_method_row_list[[length(nonbootstrap_method_row_list) + 1L]] <- build_method_row(
      dataset_key = current_key,
      source_dataset = current_source,
      analysis_variant = current_variant,
      method_name = "maller_zhou_alpha_n",
      method_role = "descriptive_tail_summary",
      statistic_name = "alpha_n",
      statistic_value = mz_object$result$alpha_n,
      p_value = mz_object$result$alpha_n,
      p_value_type = "estimated_alpha_n",
      method_flag = classify_alpha_n(mz_object$result$alpha_n),
      alpha = alpha_screening,
      tau_year = current_tau,
      mz_alpha_n = mz_object$result$alpha_n,
      mz_q_n = mz_object$result$q_n,
      convergence_flag = TRUE,
      implementation_label = "extreme_interval_followup_summary",
      note = "Maller-Zhou alpha_n and q_n are exported as descriptive tail summaries only."
    )
    
    xie_observed_object <- compute_xie_observed_object(
      df = current_data,
      dataset_key = current_key,
      source_dataset = current_source,
      analysis_variant = current_variant,
      seed_base = bootstrap_seed + i * 1000L
    )
    safe_save_rds_atomic(xie_observed_object, make_xie_observed_file(current_key))
    fit_objects_base$xie_observed[[current_key]] <- xie_observed_object
    
    receus_object <- run_receus_suite(current_data)
    fit_objects_base$receus[[current_key]] <- receus_object
    
    receus_candidate_list[[length(receus_candidate_list) + 1L]] <- receus_object$candidate_table %>%
      mutate(dataset_key = current_key, source_dataset = current_source, analysis_variant = current_variant) %>%
      select(
        dataset_key, source_dataset, analysis_variant,
        family, converged, loglik, aic,
        cure_fraction_hat, susceptible_survival_tau,
        overall_survival_tau, receus_ratio_hat,
        tau_year, receus_flag, selected_by_aic, note
      )
    
    nonbootstrap_method_row_list[[length(nonbootstrap_method_row_list) + 1L]] <- build_method_row(
      dataset_key = current_key,
      source_dataset = current_source,
      analysis_variant = current_variant,
      method_name = "receus_weibull",
      method_role = "joint_screen_sensitivity",
      statistic_name = "r_n",
      statistic_value = extract_list_scalar(receus_object$weibull_fit, "receus_ratio_hat", NA_real_),
      p_value = NA_real_,
      p_value_type = "threshold_rule",
      method_flag = classify_receus_flag(
        extract_list_scalar(receus_object$weibull_fit, "cure_fraction_hat", NA_real_),
        extract_list_scalar(receus_object$weibull_fit, "receus_ratio_hat", NA_real_)
      ),
      alpha = alpha_screening,
      model_family = "weibull",
      selected_family = "weibull",
      tau_year = current_tau,
      cure_fraction_hat = extract_list_scalar(receus_object$weibull_fit, "cure_fraction_hat", NA_real_),
      receus_ratio_hat = extract_list_scalar(receus_object$weibull_fit, "receus_ratio_hat", NA_real_),
      overall_survival_tau = extract_list_scalar(receus_object$weibull_fit, "overall_survival_tau", NA_real_),
      susceptible_survival_tau = extract_list_scalar(receus_object$weibull_fit, "susceptible_survival_tau", NA_real_),
      aic = extract_list_scalar(receus_object$weibull_fit, "aic", NA_real_),
      convergence_flag = isTRUE(extract_list_scalar(receus_object$weibull_fit, "converged", FALSE)),
      implementation_label = "mle_threshold_rule_fixed_family",
      note = "RECeUS with a fixed Weibull mixture-cure working model."
    )
    
    nonbootstrap_method_row_list[[length(nonbootstrap_method_row_list) + 1L]] <- build_method_row(
      dataset_key = current_key,
      source_dataset = current_source,
      analysis_variant = current_variant,
      method_name = "receus_aic",
      method_role = "primary_joint_appropriateness_gate",
      statistic_name = "r_n",
      statistic_value = extract_list_scalar(receus_object$best_fit, "receus_ratio_hat", NA_real_),
      p_value = NA_real_,
      p_value_type = "threshold_rule",
      method_flag = classify_receus_flag(
        extract_list_scalar(receus_object$best_fit, "cure_fraction_hat", NA_real_),
        extract_list_scalar(receus_object$best_fit, "receus_ratio_hat", NA_real_)
      ),
      alpha = alpha_screening,
      model_family = extract_list_scalar(receus_object$best_fit, "family", NA_character_),
      selected_family = receus_object$best_family,
      tau_year = current_tau,
      cure_fraction_hat = extract_list_scalar(receus_object$best_fit, "cure_fraction_hat", NA_real_),
      receus_ratio_hat = extract_list_scalar(receus_object$best_fit, "receus_ratio_hat", NA_real_),
      overall_survival_tau = extract_list_scalar(receus_object$best_fit, "overall_survival_tau", NA_real_),
      susceptible_survival_tau = extract_list_scalar(receus_object$best_fit, "susceptible_survival_tau", NA_real_),
      aic = extract_list_scalar(receus_object$best_fit, "aic", NA_real_),
      convergence_flag = isTRUE(extract_list_scalar(receus_object$best_fit, "converged", FALSE)),
      implementation_label = "mle_threshold_rule_aic_selected_family",
      note = if (!is.na(receus_object$best_family)) {
        sprintf("RECeUS-AIC using the selected family `%s`.", receus_object$best_family)
      } else {
        "RECeUS-AIC unavailable because no candidate family converged."
      }
    )
  }
  
  for (dataset_idx in seq_along(all_variant_map)) {
    current_key <- names(all_variant_map)[[dataset_idx]]
    current_map <- all_variant_map[[dataset_idx]]
    current_data <- current_map$data
    current_source <- current_map$source_dataset
    current_variant <- current_map$analysis_variant
    current_tau <- max(as.numeric(current_data$time_year))
    
    fit_objects_base$hsu_observed[[current_key]] <- list()
    
    for (family_idx in seq_along(candidate_hsu_families)) {
      family_name <- candidate_hsu_families[[family_idx]]
      log_step("Observed-only HSU setup: ", current_key, " / ", family_name)
      hsu_observed_result <- compute_hsu_supscore_from_data(
        df = current_data,
        dataset_key = current_key,
        formula_registry = ingest$formula_registry,
        working_family = family_name
      )
      hsu_observed_object <- list(
        dataset_key = current_key,
        source_dataset = current_source,
        analysis_variant = current_variant,
        working_family = family_name,
        seed_base = as.integer(bootstrap_seed + dataset_idx * 100000L + family_idx * 1000L),
        tau_year = current_tau,
        result = hsu_observed_result
      )
      safe_save_rds_atomic(hsu_observed_object, make_hsu_observed_file(current_key, family_name))
      fit_objects_base$hsu_observed[[current_key]][[family_name]] <- hsu_observed_object
    }
  }
  
  fit_objects_base$inheritance_notes[["merged__site_adjusted"]] <- "Raw-tail Maller-Zhou, Xie, and RECeUS objects are inherited conceptually from merged__site_free; only HSU is rerun with site adjustment."
  
  followup_summary_base <- compact_bind_rows(followup_summary_list) %>%
    arrange(match(dataset_key, dataset_order_full))
  
  nonbootstrap_method_results_base <- compact_bind_rows(nonbootstrap_method_row_list) %>%
    arrange(match(dataset_key, dataset_order_full), match(method_name, method_order))
  
  receus_candidate_base <- compact_bind_rows(receus_candidate_list) %>%
    arrange(match(dataset_key, dataset_order_full), match(family, candidate_receus_families))
  
  shard_registry <- build_shard_registry()
  
  stage6_run_manifest <- build_stage6_run_manifest()
  safe_write_csv_atomic(stage6_run_manifest, run_manifest_file)
  safe_write_csv_atomic(variant_registry, file.path(export_path, "stage6_variant_registry.csv"))
  write_shard_registry(shard_registry)
  
  state_bundle <- list(
    stage = "Stage 6",
    created_at = now_utc_chr(),
    session_info = utils::sessionInfo(),
    config = list(
      run_mode = run_mode,
      run_root_dir = run_root_dir,
      export_path = export_path,
      canonical_export_path = canonical_export_path,
      stage6_specification_file = stage6_specification_file,
      bootstrap_reps_xie = bootstrap_reps_xie,
      bootstrap_reps_hsu = bootstrap_reps_hsu,
      xie_shard_size = xie_shard_size,
      hsu_shard_size = hsu_shard_size,
      alpha_screening = alpha_screening,
      receus_pi_threshold = receus_pi_threshold,
      receus_r_threshold = receus_r_threshold,
      receus_pi_equivocal = receus_pi_equivocal,
      receus_r_equivocal = receus_r_equivocal,
      hsu_family_set_name = hsu_family_set_name,
      candidate_receus_families = candidate_receus_families,
      candidate_hsu_families = candidate_hsu_families,
      task_stale_after_sec = task_stale_after_sec,
      failed_task_retry_after_sec = failed_task_retry_after_sec,
      max_task_attempts = max_task_attempts,
      lock_stale_after_sec = lock_stale_after_sec
    ),
    inputs = list(
      formula_registry = ingest$formula_registry,
      common_horizons_year = ingest$common_horizons_year,
      risk_thresholds = ingest$risk_thresholds,
      source_lookup = ingest$source_lookup
    ),
    registries = list(
      variant_registry = variant_registry
    ),
    observed_outputs = list(
      followup_summary_base = followup_summary_base,
      nonbootstrap_method_results_base = nonbootstrap_method_results_base,
      receus_candidate_base = receus_candidate_base,
      fit_objects_base = fit_objects_base
    )
  )
  
  safe_save_rds_atomic(state_bundle, state_bundle_file)
  writeLines(now_utc_chr(), state_initialized_file)
  log_step("Stage 6 observed-only setup and shard registry prepared successfully.")
  write_runtime_status("state_initialized", note = "Observed-only setup completed and shard registry created.")
  invisible(TRUE)
}

initialize_if_needed <- function() {
  with_dir_lock(init_lock_dir, {
    if (isTRUE(refresh_run_state) || !file.exists(state_initialized_file) || !file.exists(state_bundle_file) || !file.exists(shard_registry_file)) {
      prepare_observed_and_registry()
    } else {
      registry_df <- requeue_stale_and_retryable_tasks(read_shard_registry())
      write_shard_registry(registry_df)
    }
  }, timeout_sec = 600, poll_sec = 0.5)
  invisible(TRUE)
}

build_hsu_method_rows_from_results <- function(dataset_key, source_dataset, analysis_variant, tau_year, family_results, familyset_result) {
  out_rows <- list()
  
  for (family_name in candidate_hsu_families) {
    obj <- family_results[[family_name]]
    adjusted_p_val <- extract_frame_cell(
      familyset_result$familywise_table[familyset_result$familywise_table$working_family == family_name, , drop = FALSE],
      "adjusted_p_value",
      NA_real_
    )
    
    out_rows[[length(out_rows) + 1L]] <- build_method_row(
      dataset_key = dataset_key,
      source_dataset = source_dataset,
      analysis_variant = analysis_variant,
      method_name = hsu_method_name(family_name),
      method_role = "presence_only_supporting_test",
      statistic_name = "sup_lr_boot",
      statistic_value = obj$statistic,
      p_value = obj$p_value,
      adjusted_p_value = adjusted_p_val,
      p_value_type = "parametric_bootstrap_approx",
      method_flag = classify_presence_p(obj$p_value),
      alpha = alpha_screening,
      model_family = family_name,
      selected_family = family_name,
      family_set_name = hsu_family_set_name,
      latency_formula_rhs = strip_formula_rhs(obj$latency_formula),
      incidence_formula_rhs = strip_formula_rhs(obj$selected_formula),
      bootstrap_reps = bootstrap_reps_hsu,
      tau_year = tau_year,
      convergence_flag = isTRUE(extract_list_scalar(obj$null_fit, "converged", FALSE)),
      implementation_label = "approximate_bootstrap_supremum_lr_family_specific",
      note = obj$note
    )
  }
  
  selected_family_final <- familyset_result$familyset$selected_family
  selected_latency_formula_rhs <- if (!is.na(selected_family_final) && nzchar(selected_family_final) && !is.null(family_results[[selected_family_final]])) {
    strip_formula_rhs(family_results[[selected_family_final]]$latency_formula)
  } else {
    NA_character_
  }
  
  out_rows[[length(out_rows) + 1L]] <- build_method_row(
    dataset_key = dataset_key,
    source_dataset = source_dataset,
    analysis_variant = analysis_variant,
    method_name = hsu_familyset_method_name,
    method_role = "presence_only_supporting_test",
    statistic_name = "familyset_sup_lr",
    statistic_value = familyset_result$familyset$statistic,
    p_value = familyset_result$familyset$p_value,
    adjusted_p_value = familyset_result$familyset$adjusted_p_value,
    p_value_type = "holm_adjusted_familywise",
    method_flag = classify_presence_p(familyset_result$familyset$p_value),
    alpha = alpha_screening,
    model_family = "familyset",
    selected_family = familyset_result$familyset$selected_family,
    family_set_name = familyset_result$familyset$family_set_name,
    latency_formula_rhs = selected_latency_formula_rhs,
    incidence_formula_rhs = strip_formula_rhs(familyset_result$familyset$selected_formula),
    bootstrap_reps = bootstrap_reps_hsu,
    tau_year = tau_year,
    convergence_flag = TRUE,
    implementation_label = "holm_adjusted_familyset_from_family_specific_hsu",
    note = familyset_result$familyset$note
  )
  
  compact_bind_rows(out_rows)
}

build_xie_method_row_from_result <- function(dataset_key, source_dataset, analysis_variant, tau_year, xie_result) {
  build_method_row(
    dataset_key = dataset_key,
    source_dataset = source_dataset,
    analysis_variant = analysis_variant,
    method_name = "xie_sufficient_followup",
    method_role = "followup_contradiction_check",
    statistic_name = "T_n",
    statistic_value = xie_result$observed$T_n,
    p_value = xie_result$p_centered,
    p_value_type = "centered_bootstrap",
    method_flag = classify_xie_flag(xie_result$p_centered),
    alpha = alpha_screening,
    bootstrap_reps = bootstrap_reps_xie,
    tau_year = tau_year,
    noncure_hat_naive = xie_result$observed$p_hat_n,
    noncure_hat_corrected = xie_result$observed$p_hat_G,
    xie_epsilon = xie_result$observed$epsilon,
    convergence_flag = TRUE,
    implementation_label = "sharded_bootstrap_extremes_test",
    note = "No-covariate Xie sufficient follow-up test with sharded centered bootstrap p-value."
  )
}

build_stage6_metadata_registry <- function(common_horizons_year, risk_thresholds) {
  tibble::tibble(
    metadata_group = c(
      "stage", "stage", "stage", "stage",
      "documents", "documents", "documents",
      "inputs", "inputs",
      "grids", "grids",
      "screening", "screening", "screening", "screening", "screening", "screening", "screening", "screening",
      "screening", "screening", "screening", "screening", "screening", "screening", "screening",
      "runtime", "runtime", "runtime", "runtime", "runtime",
      "outputs", "outputs", "outputs", "outputs"
    ),
    metadata_name = c(
      "stage_name", "stage_role", "screening_repeated_in_stage8", "run_mode",
      "canonical_common_rules", "canonical_framework", "stage6_specification_file",
      "stage1_bundle_file", "stage1_datasets_file",
      "horizon_vector", "threshold_vector",
      "alpha_screening", "xie_bootstrap_reps", "hsu_bootstrap_reps", "xie_shard_size", "hsu_shard_size",
      "receus_candidate_families", "hsu_candidate_families", "hsu_family_set_name",
      "hsu_familyset_aggregation", "decision_rule_version", "primary_gate_method", "decision_rule_presence_modifier",
      "decision_rule_xie_modifier", "decision_rule_alpha_n_role", "screening_context",
      "parallel_cores", "allow_parallel_bootstrap", "shard_registry_file", "progress_log_file", "runtime_status_file",
      "save_folder_rule", "carry_forward_join_key", "canonical_export_path", "copy_final_exports_to_canonical"
    ),
    metadata_value = c(
      "Stage 6 cure appropriateness screening",
      "Pre-cure screening and carry-forward eligibility gate only",
      "FALSE",
      run_mode,
      "2.General Model Specifications_🇬🇧ENG.md",
      "5.Model Specification Framework_🇬🇧ENG.md",
      stage6_specification_file,
      stage1_bundle_file,
      stage1_datasets_file,
      paste(common_horizons_year, collapse = "|"),
      paste(format(risk_thresholds, trim = TRUE, scientific = FALSE), collapse = "|"),
      as.character(alpha_screening),
      as.character(bootstrap_reps_xie),
      as.character(bootstrap_reps_hsu),
      as.character(xie_shard_size),
      as.character(hsu_shard_size),
      paste(candidate_receus_families, collapse = "|"),
      paste(candidate_hsu_families, collapse = "|"),
      hsu_family_set_name,
      "Holm-adjusted familywise fallback across family-specific sharded HSU working families",
      decision_rule_version,
      "RECeUS-AIC",
      "Presence-only tests can upgrade an unsupportive primary gate to equivocal but cannot create supportive status on their own",
      "Xie contradiction can downgrade a supportive primary gate to equivocal but cannot create supportive status on its own",
      "Maller-Zhou alpha_n is descriptive only and never serves as a final positive gate",
      screening_context,
      as.character(parallel_cores),
      as.character(allow_parallel_bootstrap),
      shard_registry_file,
      progress_log_file,
      runtime_status_file,
      "Write all outputs into the single user-specified run_root_dir without creating extra subfolders for canonical exports unless explicitly requested.",
      "dataset_key",
      canonical_export_path,
      as.character(copy_final_exports_to_canonical)
    )
  )
}

write_stage6_exports <- function(
    state_bundle,
    followup_sufficiency_summary,
    screening_method_results,
    receus_candidate_fits,
    screening_summary,
    carry_forward_flag_table,
    fit_objects_final,
    stage6_metadata_registry,
    stage6_shard_registry
) {
  variant_registry <- state_bundle$registries$variant_registry
  
  screening_bundle <- list(
    stage = "Stage 6",
    created_at = as.character(Sys.time()),
    session_info = utils::sessionInfo(),
    config = c(state_bundle$config, list(
      shard_registry_file = shard_registry_file
    )),
    source_lookup = state_bundle$inputs$source_lookup,
    registries = list(
      variant_registry = variant_registry,
      stage6_metadata_registry = stage6_metadata_registry,
      shard_registry = stage6_shard_registry
    ),
    outputs = list(
      followup_sufficiency_summary = followup_sufficiency_summary,
      screening_method_results = screening_method_results,
      receus_candidate_fits = receus_candidate_fits,
      screening_summary = screening_summary,
      carry_forward_flag_table = carry_forward_flag_table
    )
  )
  
  safe_write_csv_atomic(variant_registry, file.path(export_path, "stage6_variant_registry.csv"))
  safe_write_csv_atomic(stage6_metadata_registry, file.path(export_path, "stage6_metadata_registry.csv"))
  safe_write_csv_atomic(followup_sufficiency_summary, file.path(export_path, "stage6_followup_sufficiency_summary.csv"))
  safe_write_csv_atomic(screening_method_results, file.path(export_path, "stage6_screening_method_results.csv"))
  safe_write_csv_atomic(receus_candidate_fits, file.path(export_path, "stage6_receus_candidate_fits.csv"))
  safe_write_csv_atomic(screening_summary, file.path(export_path, "stage6_screening_summary.csv"))
  safe_write_csv_atomic(carry_forward_flag_table, file.path(export_path, "stage6_carry_forward_flag_table.csv"))
  safe_write_csv_atomic(stage6_shard_registry, file.path(export_path, "stage6_shard_registry.csv"))
  safe_save_rds_atomic(screening_bundle, file.path(export_path, "stage6_screening_bundle.rds"))
  
  if (isTRUE(save_fit_objects)) {
    safe_save_rds_atomic(fit_objects_final, file.path(export_path, "stage6_fitted_objects.rds"))
  }
  
  # Visual summaries
  plot_dataset_levels <- dataset_order_full
  plot_method_levels <- method_order
  plot_flag_levels <- c("supportive", "equivocal", "unsupportive")
  plot_flag_palette <- c(
    supportive = "#1B9E77",
    equivocal = "#E6AB02",
    unsupportive = "#D95F02"
  )
  plot_decision_palette <- c(
    supportive = "#1B9E77",
    equivocal = "#7570B3",
    unsupportive = "#D95F02"
  )
  plot_width_main <- 13
  plot_height_main <- 8
  plot_height_tall <- 10
  
  followup_metric_lookup <- tibble::tibble(
    metric = c(
      "tau_year",
      "largest_event_time_year",
      "plateau_length_year",
      "censor_rate_main",
      "km_tail_survival",
      "n_censored_after_last_event"
    ),
    metric_label = c(
      "Max follow-up (years)",
      "Largest event time (years)",
      "Plateau length (years)",
      "Main censoring rate",
      "KM tail survival",
      "Censored after last event"
    ),
    label_type = c("number", "number", "number", "percent", "percent", "count")
  )
  
  followup_plot_data <- dplyr::bind_rows(lapply(followup_metric_lookup$metric, function(metric_name) {
    followup_sufficiency_summary %>%
      dplyr::transmute(
        dataset_key = dataset_key,
        analysis_variant = analysis_variant,
        metric = metric_name,
        metric_value = .data[[metric_name]]
      )
  })) %>%
    dplyr::left_join(followup_metric_lookup, by = "metric") %>%
    dplyr::mutate(
      dataset_key = factor(dataset_key, levels = plot_dataset_levels),
      metric_label = factor(metric_label, levels = followup_metric_lookup$metric_label),
      metric_value_label = dplyr::case_when(
        label_type == "percent" ~ scales::label_percent(accuracy = 0.1)(metric_value),
        label_type == "count" ~ ifelse(is.na(metric_value), NA_character_, format(round(metric_value), big.mark = ",", trim = TRUE)),
        TRUE ~ ifelse(is.na(metric_value), NA_character_, format(round(metric_value, 2), nsmall = 2, trim = TRUE))
      )
    )
  
  followup_plot <- ggplot2::ggplot(followup_plot_data, ggplot2::aes(x = metric_value, y = dataset_key)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = metric_value, yend = dataset_key), linewidth = 0.45, colour = "#BDBDBD", na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(colour = analysis_variant), size = 2.8, na.rm = TRUE) +
    ggplot2::geom_text(
      ggplot2::aes(label = metric_value_label),
      hjust = -0.10,
      size = 3.1,
      na.rm = TRUE
    ) +
    ggplot2::facet_wrap(~metric_label, scales = "free_x", ncol = 2) +
    ggplot2::scale_colour_manual(
      values = c(single_cohort = "#1F78B4", site_free = "#33A02C", site_adjusted = "#E31A1C"),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Stage 6 follow-up sufficiency summary",
      subtitle = "Observed raw-tail summaries exported from the sharded Stage 6 reducer",
      x = NULL,
      y = NULL,
      colour = "Analysis variant"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  screening_method_plot_data <- screening_method_results %>%
    dplyr::mutate(
      dataset_key = factor(dataset_key, levels = plot_dataset_levels),
      method_name = factor(method_name, levels = plot_method_levels),
      method_label = factor(stringr::str_replace_all(as.character(method_name), "_", "\n"), levels = stringr::str_replace_all(plot_method_levels, "_", "\n")),
      method_flag = factor(method_flag, levels = plot_flag_levels),
      plot_label = dplyr::case_when(
        !is.na(p_value) ~ paste0("p=", formatC(p_value, format = "f", digits = 3)),
        !is.na(receus_ratio_hat) & !is.na(cure_fraction_hat) ~ paste0("r=", formatC(receus_ratio_hat, format = "f", digits = 2), "\ncf=", formatC(cure_fraction_hat, format = "f", digits = 2)),
        !is.na(statistic_value) ~ formatC(statistic_value, format = "f", digits = 2),
        TRUE ~ ""
      )
    )
  
  screening_method_plot <- ggplot2::ggplot(screening_method_plot_data, ggplot2::aes(x = method_label, y = dataset_key, fill = method_flag)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = plot_label), size = 2.6, lineheight = 0.95, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = plot_flag_palette, drop = FALSE, na.value = "#D9D9D9") +
    ggplot2::labs(
      title = "Stage 6 screening method flags",
      subtitle = "RECeUS-AIC primary gate, family-specific HSU rows, and HSU family-set row",
      x = NULL,
      y = NULL,
      fill = "Method flag"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 8)
    )
  
  receus_plot_data <- receus_candidate_fits %>%
    dplyr::mutate(
      dataset_key = factor(dataset_key, levels = plot_dataset_levels),
      family = factor(family, levels = candidate_receus_families),
      selected_status = dplyr::if_else(selected_by_aic, "AIC selected", "Candidate"),
      selected_status = factor(selected_status, levels = c("Candidate", "AIC selected"))
    )
  
  receus_plot <- ggplot2::ggplot(receus_plot_data, ggplot2::aes(x = receus_ratio_hat, y = cure_fraction_hat)) +
    ggplot2::geom_vline(xintercept = receus_r_threshold, linetype = "dashed", colour = "#636363") +
    ggplot2::geom_hline(yintercept = receus_pi_threshold, linetype = "dashed", colour = "#636363") +
    ggplot2::geom_point(ggplot2::aes(shape = family, colour = selected_status), size = 3.0, alpha = 0.9, na.rm = TRUE) +
    ggplot2::facet_wrap(~dataset_key, ncol = 2) +
    ggplot2::scale_colour_manual(values = c("Candidate" = "#6A3D9A", "AIC selected" = "#E31A1C"), drop = FALSE) +
    ggplot2::labs(
      title = "Stage 6 RECeUS candidate-family map",
      subtitle = "RECeUS-AIC remains the primary eligibility gate",
      x = "RECeUS ratio estimate",
      y = "Estimated cure fraction",
      colour = "Selection status",
      shape = "Latency family"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  carry_forward_base <- carry_forward_flag_table %>%
    dplyr::transmute(
      dataset_key = factor(dataset_key, levels = plot_dataset_levels),
      cure_model_eligibility_flag,
      presence_evidence_label,
      followup_evidence_label,
      presence_component = dplyr::case_when(
        dplyr::coalesce(as.logical(presence_support_flag), FALSE) ~ "supportive",
        dplyr::coalesce(as.logical(presence_unsupport_flag), FALSE) ~ "unsupportive",
        TRUE ~ "equivocal"
      ),
      followup_component = dplyr::case_when(
        dplyr::coalesce(as.logical(followup_support_flag), FALSE) ~ "supportive",
        dplyr::coalesce(as.logical(followup_unsupport_flag), FALSE) ~ "unsupportive",
        TRUE ~ "equivocal"
      )
    )
  
  carry_forward_plot_data <- dplyr::bind_rows(
    carry_forward_base %>%
      dplyr::transmute(
        dataset_key = dataset_key,
        component = "Presence evidence",
        decision_flag = presence_component,
        component_label = presence_evidence_label
      ),
    carry_forward_base %>%
      dplyr::transmute(
        dataset_key = dataset_key,
        component = "Primary gate / Xie check",
        decision_flag = followup_component,
        component_label = followup_evidence_label
      ),
    carry_forward_base %>%
      dplyr::transmute(
        dataset_key = dataset_key,
        component = "Eligibility flag",
        decision_flag = cure_model_eligibility_flag,
        component_label = cure_model_eligibility_flag
      )
  ) %>%
    dplyr::mutate(
      component = factor(
        component,
        levels = c("Presence evidence", "Primary gate / Xie check", "Eligibility flag")
      ),
      decision_flag = factor(decision_flag, levels = c("supportive", "equivocal", "unsupportive")),
      component_label = as.character(component_label)
    )
  
  carry_forward_plot <- ggplot2::ggplot(carry_forward_plot_data, ggplot2::aes(x = component, y = dataset_key, fill = decision_flag)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = stringr::str_replace_all(component_label, "_", "\n")), size = 2.8, lineheight = 0.92, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = plot_decision_palette, drop = FALSE, na.value = "#D9D9D9") +
    ggplot2::labs(
      title = "Stage 6 carry-forward eligibility map",
      subtitle = "Decision rule: RECeUS-AIC primary gate, presence-only modifier, Xie contradiction downgrade",
      x = NULL,
      y = NULL,
      fill = "Decision flag"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank()
    )
  
  safe_save_plot(followup_plot, file.path(export_path, "stage6_followup_sufficiency_summary.png"), width = plot_width_main, height = plot_height_tall)
  safe_save_plot(screening_method_plot, file.path(export_path, "stage6_screening_method_results.png"), width = plot_width_main, height = plot_height_main)
  safe_save_plot(receus_plot, file.path(export_path, "stage6_receus_candidate_fits.png"), width = plot_width_main, height = plot_height_main)
  safe_save_plot(carry_forward_plot, file.path(export_path, "stage6_carry_forward_flag_table.png"), width = plot_width_main, height = plot_height_main)
  
  export_file_names <- c(
    "stage6_variant_registry.csv",
    "stage6_metadata_registry.csv",
    "stage6_followup_sufficiency_summary.csv",
    "stage6_screening_method_results.csv",
    "stage6_receus_candidate_fits.csv",
    "stage6_screening_summary.csv",
    "stage6_carry_forward_flag_table.csv",
    "stage6_shard_registry.csv",
    "stage6_screening_bundle.rds",
    "stage6_followup_sufficiency_summary.png",
    "stage6_screening_method_results.png",
    "stage6_receus_candidate_fits.png",
    "stage6_carry_forward_flag_table.png",
    "stage6_run_manifest.csv",
    "stage6_progress_log.txt",
    "stage6_runtime_status.csv"
  )
  
  export_object_names <- c(
    "variant_registry",
    "stage6_metadata_registry",
    "followup_sufficiency_summary",
    "screening_method_results",
    "receus_candidate_fits",
    "screening_summary",
    "carry_forward_flag_table",
    "stage6_shard_registry",
    "screening_bundle",
    "followup_plot",
    "screening_method_plot",
    "receus_plot",
    "carry_forward_plot",
    "stage6_run_manifest",
    "progress_log",
    "runtime_status"
  )
  
  export_descriptions <- c(
    "Stage 6 structural variant registry keyed by dataset_key.",
    "Stage 6 metadata and carry-forward rules.",
    "Dataset-level follow-up sufficiency summary including extreme-tail quantities.",
    "Long-format Stage 6 screening method results including family-specific HSU rows and HSU family-set row.",
    "RECeUS and RECeUS-AIC candidate-family fit summary table.",
    "Wide Stage 6 screening summary by dataset_key.",
    "Dataset-key-based carry-forward eligibility flag table for later stages.",
    "Shard registry tracking pending/running/completed shard tasks.",
    "Reusable Stage 6 bundle with summary tables and inherited backbone settings.",
    "Visualization exported from stage6_followup_sufficiency_summary.csv.",
    "Visualization exported from stage6_screening_method_results.csv.",
    "Visualization exported from stage6_receus_candidate_fits.csv.",
    "Visualization exported from stage6_carry_forward_flag_table.csv.",
    "Run manifest describing the sharded Stage 6 runtime.",
    "Console and file progress log written during the run.",
    "Runtime status log written during the run."
  )
  
  if (isTRUE(save_fit_objects)) {
    export_file_names <- c(export_file_names, "stage6_fitted_objects.rds")
    export_object_names <- c(export_object_names, "fit_objects")
    export_descriptions <- c(export_descriptions, "Observed objects, reduced bootstrap distributions, and final Stage 6 fitted-object bundle.")
  }
  
  stage6_export_manifest <- tibble::tibble(
    file_name = export_file_names,
    object_name = export_object_names,
    description = export_descriptions,
    file_path = file.path(export_path, export_file_names)
  )
  
  safe_write_csv_atomic(stage6_export_manifest, file.path(export_path, "stage6_export_manifest.csv"))
  
  copy_final_outputs_to_canonical(c(stage6_export_manifest$file_name, "stage6_export_manifest.csv"))
  
  invisible(TRUE)
}

reduce_and_export_stage6 <- function() {
  state_bundle <- readRDS(state_bundle_file)
  registry_df <- read_shard_registry()
  
  if (any(registry_df$status %in% c("pending", "running"))) {
    return(FALSE)
  }
  
  if (any(registry_df$status == "failed")) {
    log_step("Reduction blocked because failed shard(s) remain in the registry.", level = "WARN")
    return(FALSE)
  }
  
  log_step("All shards complete. Starting Stage 6 reducer.")
  write_runtime_status("reducing", note = "All shards complete; Stage 6 reducer started.")
  
  variant_registry <- state_bundle$registries$variant_registry
  followup_sufficiency_summary <- state_bundle$observed_outputs$followup_summary_base
  nonbootstrap_method_results_base <- state_bundle$observed_outputs$nonbootstrap_method_results_base
  receus_candidate_fits <- state_bundle$observed_outputs$receus_candidate_base
  fit_objects_final <- state_bundle$observed_outputs$fit_objects_base
  
  base_variant_map <- build_base_variant_map(readRDS(prepared_datasets_file))
  all_variant_map <- build_dataset_map_all(readRDS(prepared_datasets_file))
  
  xie_method_row_list <- list()
  xie_results_final <- list()
  
  for (dataset_key in names(base_variant_map)) {
    observed_object <- readRDS(make_xie_observed_file(dataset_key))
    bootstrap_stats <- collect_bootstrap_stats_from_registry(
      registry_df = registry_df,
      task_type = "xie",
      dataset_key = dataset_key,
      total_reps = bootstrap_reps_xie
    )
    xie_result <- finalize_xie_from_bootstrap_stats(observed_object, bootstrap_stats)
    xie_results_final[[dataset_key]] <- xie_result
    
    xie_method_row_list[[length(xie_method_row_list) + 1L]] <- build_xie_method_row_from_result(
      dataset_key = dataset_key,
      source_dataset = observed_object$source_dataset,
      analysis_variant = observed_object$analysis_variant,
      tau_year = observed_object$tau_year,
      xie_result = xie_result
    )
  }
  
  xie_method_results_base <- compact_bind_rows(xie_method_row_list)
  
  hsu_method_row_list <- list()
  hsu_results_final <- list()
  
  for (dataset_key in names(all_variant_map)) {
    current_map <- all_variant_map[[dataset_key]]
    current_tau <- max(as.numeric(current_map$data$time_year))
    family_results <- list()
    
    for (family_name in candidate_hsu_families) {
      observed_object <- readRDS(make_hsu_observed_file(dataset_key, family_name))
      bootstrap_stats <- collect_bootstrap_stats_from_registry(
        registry_df = registry_df,
        task_type = "hsu",
        dataset_key = dataset_key,
        working_family = family_name,
        total_reps = bootstrap_reps_hsu
      )
      family_results[[family_name]] <- finalize_hsu_family_from_bootstrap_stats(observed_object, bootstrap_stats)
    }
    
    familyset_result <- reduce_hsu_familyset_from_family_results(family_results)
    
    hsu_results_final[[dataset_key]] <- list(
      family_results = family_results,
      familyset = familyset_result
    )
    
    hsu_method_row_list[[length(hsu_method_row_list) + 1L]] <- build_hsu_method_rows_from_results(
      dataset_key = dataset_key,
      source_dataset = current_map$source_dataset,
      analysis_variant = current_map$analysis_variant,
      tau_year = current_tau,
      family_results = family_results,
      familyset_result = familyset_result
    )
  }
  
  hsu_method_results_all <- compact_bind_rows(hsu_method_row_list)
  
  screening_method_results_base <- bind_rows(
    nonbootstrap_method_results_base,
    xie_method_results_base,
    hsu_method_results_all %>% filter(dataset_key != "merged__site_adjusted")
  ) %>%
    arrange(match(dataset_key, dataset_order_full), match(method_name, method_order))
  
  merged_followup_copy <- followup_sufficiency_summary %>%
    filter(dataset_key == "merged__site_free") %>%
    mutate(
      dataset_key = "merged__site_adjusted",
      analysis_variant = "site_adjusted",
      summary_note = append_note(
        summary_note,
        "Raw follow-up, Maller-Zhou, Xie, and RECeUS tail metrics are inherited from the merged site-free data because the underlying merged cohort is unchanged; only the HSU working-model family set is rerun with site adjustment."
      )
    )
  
  followup_sufficiency_summary <- bind_rows(followup_sufficiency_summary, merged_followup_copy) %>%
    arrange(match(dataset_key, dataset_order_full))
  
  methods_inherited_to_site_adjusted <- c("maller_zhou_dn", "maller_zhou_alpha_n", "xie_sufficient_followup", "receus_weibull", "receus_aic")
  
  merged_inherited_methods <- screening_method_results_base %>%
    filter(dataset_key == "merged__site_free", method_name %in% methods_inherited_to_site_adjusted) %>%
    mutate(
      dataset_key = "merged__site_adjusted",
      analysis_variant = "site_adjusted",
      note = append_note(note, "Raw-tail metric inherited from the merged site-free branch for the merged site-adjusted structural variant.")
    )
  
  screening_method_results <- bind_rows(
    screening_method_results_base,
    merged_inherited_methods,
    hsu_method_results_all %>% filter(dataset_key == "merged__site_adjusted")
  ) %>%
    arrange(match(dataset_key, dataset_order_full), match(method_name, method_order))
  
  merged_inherited_receus <- receus_candidate_fits %>%
    filter(dataset_key == "merged__site_free") %>%
    mutate(
      dataset_key = "merged__site_adjusted",
      analysis_variant = "site_adjusted",
      note = append_note(note, "Candidate RECeUS fits inherited from merged site-free because RECeUS is evaluated on the same merged raw tail.")
    )
  
  receus_candidate_fits <- bind_rows(receus_candidate_fits, merged_inherited_receus) %>%
    arrange(match(dataset_key, dataset_order_full), match(family, candidate_receus_families))
  
  fit_objects_final$xie <- xie_results_final
  fit_objects_final$hsu <- hsu_results_final
  
  decision_rows <- bind_rows(lapply(seq_len(nrow(variant_registry)), function(i) {
    current_variant <- variant_registry[i, , drop = FALSE]
    current_key <- current_variant$dataset_key[[1]]
    
    derive_dataset_decision_row(
      dataset_key = current_key,
      source_dataset = current_variant$source_dataset[[1]],
      analysis_variant = current_variant$analysis_variant[[1]],
      method_df = screening_method_results[screening_method_results$dataset_key == current_key, , drop = FALSE],
      common_horizons_year = state_bundle$inputs$common_horizons_year,
      risk_thresholds = state_bundle$inputs$risk_thresholds,
      screening_note = "Bayesian Stage 8 must carry this Stage 6 cure-model eligibility result forward and must not repeat cure-appropriateness screening."
    )
  }))
  
  carry_forward_flag_table <- variant_registry %>%
    left_join(decision_rows, by = c("dataset_key", "source_dataset", "analysis_variant", "screening_context")) %>%
    mutate(screening_note = append_note(screening_note, screening_variant_note)) %>%
    select(
      dataset_key, source_dataset, source_description, analysis_variant, screening_context,
      tail_metric_source, hsu_formula_branch,
      hsu_latency_formula_rhs, hsu_incidence_formula_rhs_set, hsu_working_family_set,
      decision_rule_version, primary_gate_method,
      primary_gate_flag, presence_modifier_flag, followup_contradiction_flag,
      descriptive_tail_summary_flag, presence_modifier_label, followup_contradiction_label,
      modifier_evidence_label,
      cure_model_eligibility_flag, final_decision_flag,
      presence_support_flag, followup_support_flag,
      presence_unsupport_flag, followup_unsupport_flag,
      presence_evidence_label, followup_evidence_label,
      maller_zhou_dn_p_value, hsu_familyset_p_value,
      xie_centered_bootstrap_p_value,
      receus_weibull_flag, receus_aic_flag,
      supporting_methods, contradicting_methods,
      common_horizon_vector, common_threshold_vector,
      carry_forward_stage7, carry_forward_stage8,
      screening_repeated_in_stage8,
      screening_note
    ) %>%
    arrange(match(dataset_key, dataset_order_full))
  
  screening_summary <- bind_rows(lapply(dataset_order_full, function(current_key) {
    build_screening_summary_row(
      dataset_key = current_key,
      followup_summary = followup_sufficiency_summary,
      method_results = screening_method_results,
      carry_forward_table = carry_forward_flag_table
    )
  })) %>%
    arrange(match(dataset_key, dataset_order_full))
  
  stage6_metadata_registry <- build_stage6_metadata_registry(
    common_horizons_year = state_bundle$inputs$common_horizons_year,
    risk_thresholds = state_bundle$inputs$risk_thresholds
  )
  
  write_stage6_exports(
    state_bundle = state_bundle,
    followup_sufficiency_summary = followup_sufficiency_summary,
    screening_method_results = screening_method_results,
    receus_candidate_fits = receus_candidate_fits,
    screening_summary = screening_summary,
    carry_forward_flag_table = carry_forward_flag_table,
    fit_objects_final = fit_objects_final,
    stage6_metadata_registry = stage6_metadata_registry,
    stage6_shard_registry = registry_df
  )
  
  writeLines(as.character(Sys.time()), reduction_complete_file)
  write_runtime_status("completed", note = "Stage 6 sharded run completed successfully.")
  log_step("Stage 6 sharded run completed successfully.")
  TRUE
}

maybe_reduce_stage6 <- function() {
  if (file.exists(reduction_complete_file)) {
    return(TRUE)
  }
  
  reduced <- tryCatch(
    {
      with_dir_lock(reduce_lock_dir, {
        if (file.exists(reduction_complete_file)) {
          return(TRUE)
        }
        reduce_and_export_stage6()
      }, timeout_sec = worker_nonblocking_reduce_lock_timeout_sec, poll_sec = 0.25)
    },
    error = function(e) {
      FALSE
    }
  )
  
  isTRUE(reduced)
}

run_worker_loop <- function() {
  idle_cycles <- 0L
  touch_worker_heartbeat(stage = "worker_started", note = "Worker loop entered.")
  
  repeat {
    if (file.exists(reduction_complete_file)) {
      log_step("Reduction already complete; worker exiting.")
      touch_worker_heartbeat(stage = "worker_exiting", note = "Reduction already complete.")
      break
    }
    
    touch_worker_heartbeat(stage = "idle", note = paste("Idle cycle", idle_cycles))
    
    task_row <- tryCatch(claim_next_task(), error = function(e) {
      log_step("Task claim error: ", conditionMessage(e), level = "WARN")
      NULL
    })
    
    if (!is.null(task_row) && nrow(task_row) == 1L) {
      idle_cycles <- 0L
      run_one_claimed_task(task_row)
      next
    }
    
    reduced <- maybe_reduce_stage6()
    if (isTRUE(reduced) || file.exists(reduction_complete_file)) {
      touch_worker_heartbeat(stage = "worker_exiting", note = "Reduction completed.")
      break
    }
    
    idle_cycles <- idle_cycles + 1L
    if (idle_cycles >= max_worker_idle_cycles) {
      log_step(
        "Worker idle limit reached without reduction completion. Exiting cleanly after ",
        idle_cycles, " idle cycles."
      )
      touch_worker_heartbeat(stage = "worker_exiting", note = "Idle limit reached.")
      break
    }
    
    if (file.exists(shard_registry_file)) {
      registry_df <- tryCatch(read_shard_registry(), error = function(e) tibble::tibble())
      if (nrow(registry_df) > 0L && any(registry_df$status == "failed")) {
        failed_attempts <- as_integer_or_na(registry_df$attempt_count)
        failed_attempts[is.na(failed_attempts)] <- 0L
        unretryable_failed <- registry_df$status == "failed" & failed_attempts >= max_task_attempts
        if (any(unretryable_failed, na.rm = TRUE)) {
          log_step("Unretryable failed shard(s) detected in registry. Exiting worker so failures can be inspected.", level = "WARN")
          touch_worker_heartbeat(stage = "worker_exiting", note = "Unretryable failed shard detected.")
          break
        }
      }
    }
    
    Sys.sleep(worker_poll_seconds)
  }
  
  invisible(TRUE)
}

# 🔴 Main: initialize once, then let every worker auto-claim shard tasks ===============================
initialize_if_needed()
run_worker_loop()