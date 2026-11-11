## =============================================================================
##  dev/verify_metabolome_v2.R
##  ---------------------------------------------------------------------------
##  PURPOSE     One-off check after the fecal metabolite table was replaced by
##              the laboratory's corrected report (November 2026): re-run, in
##              fresh R processes, every script whose results depend on
##              input/metabolome/fecal_metabolites.csv - except the MOFA
##              scripts 20-23, which need a retrained model - and keep the
##              previous results next to the new ones for comparison.
##
##  USAGE       from the package root (the folder with LARC_TNT_response.Rproj):
##                  source("dev/verify_metabolome_v2.R")
##              or  Rscript dev/verify_metabolome_v2.R
##              (from another working directory, source() it with its full path)
##
##  WHAT IT DOES
##    1. copies output/tables, output/logs and output/figures to output_v1/
##       (only if output_v1/ does not exist yet), so the results of the first
##       laboratory report stay available;
##    2. runs scripts 07, 08, 09a, 09b, 09c, 10, 14, 15, 16, 17, 18, 19 with
##       their default switches (logs in output/logs/);
##    3. optionally refits the metabolite random forest only (all other assays
##       are replayed) from a temporary copy of script 09b with
##       refit_assays <- "Metabolite", then redraws Supplementary Figure 4.
##       Nothing in input/ is changed (freeze_refit stays FALSE).
##
##  OUTPUT      output/logs/_verify_metabolome_v2_status.csv  (script, status, seconds)
## =============================================================================

local({      # nothing is left behind in the calling R session

  refit_metabolite_rf <- TRUE     # <-- step 3 (needs randomForest, rsample, fastshap; skipped if not installed)

  ## When the working directory is not the package root, find it from this file's own path
  if (!file.exists("LARC_TNT_response.Rproj")) {
    here <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
    if (is.null(here)) {
      arg <- commandArgs(FALSE); arg <- sub("^--file=", "", arg[grepl("^--file=", arg)])
      if (length(arg) > 0) here <- arg[1]
    }
    if (!is.null(here)) setwd(dirname(dirname(normalizePath(here))))
  }
  if (!file.exists("LARC_TNT_response.Rproj")) stop("Run from the package root (the folder with LARC_TNT_response.Rproj).")
  message("Package root: ", getwd())

  ## ---- 1. keep the previous results ---------------------------------------------------
  if (dir.exists("output") && !dir.exists("output_v1")) {
    dir.create("output_v1")
    for (d in c("tables", "logs", "figures")) {
      if (dir.exists(file.path("output", d))) file.copy(file.path("output", d), "output_v1", recursive = TRUE)
    }
    message("Previous results copied to output_v1/")
  }
  dir.create("output/logs", recursive = TRUE, showWarnings = FALSE)

  ## ---- 2. run the metabolite-dependent scripts ----------------------------------------
  ids <- c("07", "08", "09a", "09b", "09c", "10", "14", "15", "16", "17", "18", "19")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  scripts <- sort(list.files("R", pattern = "^[0-9].*\\.R$"))
  scripts <- scripts[sub("_.*$", "", scripts) %in% ids]

  run_script <- function(path, log_file) {
    t0 <- Sys.time()
    code <- system2(rscript, c("--encoding=UTF-8", shQuote(path)), stdout = log_file, stderr = log_file)
    list(status = if (code == 0) "ok" else "error", seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))
  }

  status <- data.frame(script = scripts, status = NA_character_, seconds = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(scripts)) {
    s <- scripts[i]
    message(sprintf("[%2d/%d] RUN   %s", i, length(scripts), s), appendLF = FALSE)
    log_file <- file.path("output/logs", sub("\\.R$", ".log", s))
    res <- run_script(file.path("R", s), log_file)
    status$status[i] <- res$status; status$seconds[i] <- res$seconds
    message(sprintf("  %s (%.0f s)", if (res$status == "ok") "done" else paste("FAILED - see", log_file), res$seconds))
  }

  ## ---- 3. optional: refit the metabolite random forest ----------------------------------
  rf_packages <- c("randomForest", "rsample", "fastshap")
  if (refit_metabolite_rf) {
    missing_pkg <- rf_packages[!vapply(rf_packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkg) > 0) {
      message("Metabolite random-forest refit skipped (not installed: ", paste(missing_pkg, collapse = ", "), ")")
      status <- rbind(status, data.frame(script = "09b refit Metabolite", status = "skipped", seconds = NA_real_))
    } else {
      src <- readLines("R/09b_SFig4_Random_forest_models.R", encoding = "UTF-8")
      hit <- grep("^refit_assays <- character\\(0\\)", src)
      stopifnot(length(hit) == 1)
      src[hit] <- 'refit_assays <- "Metabolite"'
      tmp <- "R/_verify_09b_refit_metabolite.R"
      writeLines(src, tmp, useBytes = TRUE)
      message("[refit] RUN   09b with refit_assays = \"Metabolite\"", appendLF = FALSE)
      res <- run_script(tmp, "output/logs/_verify_09b_refit_metabolite.log")
      unlink(tmp)
      message(sprintf("  %s (%.0f s)", if (res$status == "ok") "done" else "FAILED - see output/logs/_verify_09b_refit_metabolite.log", res$seconds))
      status <- rbind(status, data.frame(script = "09b refit Metabolite", status = res$status, seconds = res$seconds))
      if (res$status == "ok") {
        s <- scripts[sub("_.*$", "", scripts) == "09c"]
        res <- run_script(file.path("R", s), "output/logs/_verify_09c_after_refit.log")
        status <- rbind(status, data.frame(script = "09c after refit", status = res$status, seconds = res$seconds))
      }
    }
  }

  write.csv(status, "output/logs/_verify_metabolome_v2_status.csv", row.names = FALSE)
  print(status, row.names = FALSE)
  message("Tables: output/tables/   logs: output/logs/   previous results: output_v1/")

})
