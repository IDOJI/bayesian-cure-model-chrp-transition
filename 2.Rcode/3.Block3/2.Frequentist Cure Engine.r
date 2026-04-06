# Configure: Block 3 frequentist cure wrapper ============================
find_repo_root <- function(start_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    has_repo_markers <- dir.exists(file.path(current_dir, ".git")) ||
      (
        dir.exists(file.path(current_dir, "0.Data")) &&
          dir.exists(file.path(current_dir, "2.Rcode")) &&
          dir.exists(file.path(current_dir, "3.Results files"))
      )

    if (has_repo_markers) {
      return(current_dir)
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }

  stop("Could not locate the repository root.", call. = FALSE)
}

command_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", command_args, value = TRUE)
search_start_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}

repo_root <- find_repo_root(search_start_dir)
shared_analysis_dir <- file.path(
  repo_root,
  "2.Rcode",
  paste0("1", ".Block1"),
  "Main"
)
engine_script <- file.path(
  shared_analysis_dir,
  paste0("3", ".Frequentist Mixture-Cure Surival Analysis (Lognormal).r")
)

Sys.setenv(
  STAGE7_SIMPLE_DATA_FILE = Sys.getenv("BLOCK3_FREQUENTIST_DATA_FILE", unset = Sys.getenv("STAGE7_SIMPLE_DATA_FILE", unset = "")),
  STAGE7_SIMPLE_EXPORT_PATH = Sys.getenv("BLOCK3_FREQUENTIST_EXPORT_PATH", unset = Sys.getenv("STAGE7_SIMPLE_EXPORT_PATH", unset = "")),
  STAGE7_SIMPLE_PNU_SITE_LABEL = Sys.getenv("BLOCK3_PNU_SITE_LABEL", unset = Sys.getenv("STAGE7_SIMPLE_PNU_SITE_LABEL", unset = "PNU")),
  STAGE7_SIMPLE_SNU_SITE_LABEL = Sys.getenv("BLOCK3_SNU_SITE_LABEL", unset = Sys.getenv("STAGE7_SIMPLE_SNU_SITE_LABEL", unset = "SNU")),
  STAGE7_SIMPLE_REUSE_EXISTING_FIT_RDS = Sys.getenv("BLOCK3_FREQUENTIST_REUSE_EXISTING_FIT_RDS", unset = Sys.getenv("STAGE7_SIMPLE_REUSE_EXISTING_FIT_RDS", unset = "TRUE")),
  STAGE7_SIMPLE_OUTPUT_PREFIX = Sys.getenv("BLOCK3_FREQUENTIST_OUTPUT_PREFIX", unset = "block3_frequentist_cure"),
  STAGE7_SIMPLE_ANALYSIS_LABEL = Sys.getenv("BLOCK3_FREQUENTIST_ANALYSIS_LABEL", unset = "Block 3 frequentist mixture cure")
)

source(engine_script, local = new.env(parent = globalenv()))
