# Configure: Block 3 Bayesian cure wrapper ===============================
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
  paste0("4", ".Bayesian Mixture Cure Model.r")
)

Sys.setenv(
  STAGE8A_SIMPLE_DATA_FILE = Sys.getenv("BLOCK3_BAYESIAN_DATA_FILE", unset = Sys.getenv("STAGE8A_SIMPLE_DATA_FILE", unset = "")),
  STAGE8A_SIMPLE_EXPORT_PATH = Sys.getenv("BLOCK3_BAYESIAN_EXPORT_PATH", unset = Sys.getenv("STAGE8A_SIMPLE_EXPORT_PATH", unset = "")),
  STAGE8A_SIMPLE_PNU_SITE_LABEL = Sys.getenv("BLOCK3_PNU_SITE_LABEL", unset = Sys.getenv("STAGE8A_SIMPLE_PNU_SITE_LABEL", unset = "PNU")),
  STAGE8A_SIMPLE_SNU_SITE_LABEL = Sys.getenv("BLOCK3_SNU_SITE_LABEL", unset = Sys.getenv("STAGE8A_SIMPLE_SNU_SITE_LABEL", unset = "SNU")),
  STAGE8A_SIMPLE_REUSE_EXISTING_FIT_RDS = Sys.getenv("BLOCK3_BAYESIAN_REUSE_EXISTING_FIT_RDS", unset = Sys.getenv("STAGE8A_SIMPLE_REUSE_EXISTING_FIT_RDS", unset = "TRUE")),
  STAGE8A_SIMPLE_STAN_CHAINS = Sys.getenv("BLOCK3_BAYESIAN_STAN_CHAINS", unset = Sys.getenv("STAGE8A_SIMPLE_STAN_CHAINS", unset = "2")),
  STAGE8A_SIMPLE_STAN_ITER = Sys.getenv("BLOCK3_BAYESIAN_STAN_ITER", unset = Sys.getenv("STAGE8A_SIMPLE_STAN_ITER", unset = "1000")),
  STAGE8A_SIMPLE_STAN_WARMUP = Sys.getenv("BLOCK3_BAYESIAN_STAN_WARMUP", unset = Sys.getenv("STAGE8A_SIMPLE_STAN_WARMUP", unset = "500")),
  STAGE8A_SIMPLE_STAN_THIN = Sys.getenv("BLOCK3_BAYESIAN_STAN_THIN", unset = Sys.getenv("STAGE8A_SIMPLE_STAN_THIN", unset = "1")),
  STAGE8A_SIMPLE_STAN_ADAPT_DELTA = Sys.getenv("BLOCK3_BAYESIAN_STAN_ADAPT_DELTA", unset = Sys.getenv("STAGE8A_SIMPLE_STAN_ADAPT_DELTA", unset = "0.95")),
  STAGE8A_SIMPLE_STAN_MAX_TREEDEPTH = Sys.getenv("BLOCK3_BAYESIAN_STAN_MAX_TREEDEPTH", unset = Sys.getenv("STAGE8A_SIMPLE_STAN_MAX_TREEDEPTH", unset = "12")),
  STAGE8A_SIMPLE_POSTERIOR_PRED_DRAWS = Sys.getenv("BLOCK3_BAYESIAN_POSTERIOR_PRED_DRAWS", unset = Sys.getenv("STAGE8A_SIMPLE_POSTERIOR_PRED_DRAWS", unset = "200")),
  STAGE8A_SIMPLE_OUTPUT_PREFIX = Sys.getenv("BLOCK3_BAYESIAN_OUTPUT_PREFIX", unset = "block3_bayesian_cure"),
  STAGE8A_SIMPLE_ANALYSIS_LABEL = Sys.getenv("BLOCK3_BAYESIAN_ANALYSIS_LABEL", unset = "Block 3 Bayesian cure")
)

source(engine_script, local = new.env(parent = globalenv()))
