## =============================================================================
##  run_all.R  --  run every analysis script in order
##  ---------------------------------------------------------------------------
##  Usage (from the package root, e.g. after opening LARC_TNT_response.Rproj):
##
##      source("run_all.R")
##
##  Every script is run in a FRESH R process (Rscript), exactly as the
##  verification runs were done, with its console output saved to
##  output/logs/<script>.log.  A fresh process per script avoids package
##  masking between scripts and keeps memory use low.  Scripts that need the
##  tumour count matrix are skipped with a message if that file has been
##  removed, so the rest of the analysis still completes.
##
##  Options
##      stop_on_error   TRUE  = abort at the first failing script
##      scripts_to_run  NULL  = all; or a vector of prefixes, e.g. c("01", "07")
##
##  Runtime (R 4.3, laptop): about 15 minutes with the default replay
##  switches; several hours if the refit switches in scripts 09b / 11 / 12 /
##  20 are turned on.
## =============================================================================

stop_on_error  <- FALSE
scripts_to_run <- NULL

if (!file.exists("LARC_TNT_response.Rproj")) stop("Run from the package root (the folder with LARC_TNT_response.Rproj).")
dir.create("output/logs", recursive = TRUE, showWarnings = FALSE)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
scripts <- sort(list.files("R", pattern = "^[0-9].*\\.R$"))
if (!is.null(scripts_to_run)) scripts <- scripts[sub("_.*$", "", scripts) %in% scripts_to_run]

## scripts that cannot run without the tumour count matrix (public, but large;
## skipped if it has been removed).  The clinical table is optional: 04b replays
## frozen results and 16 / 22 / 23 use the public covariate file.
needs <- list(
  "input/host_rnaseq/tumor_rna_counts.csv" = c("10", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23")
)

status <- data.frame(script = scripts, status = NA_character_, seconds = NA_real_, stringsAsFactors = FALSE)
for (i in seq_along(scripts)) {
  s <- scripts[i]; id <- sub("_.*$", "", s)
  missing <- names(needs)[vapply(needs, function(v) id %in% v, logical(1)) & !file.exists(names(needs))]
  if (length(missing)) {
    message(sprintf("[%2d/%d] SKIP  %s  (input file not found: %s)", i, length(scripts), s, paste(missing, collapse = ", ")))
    status$status[i] <- "skipped"; next
  }
  message(sprintf("[%2d/%d] RUN   %s", i, length(scripts), s), appendLF = FALSE)
  log_file <- file.path("output/logs", sub("\\.R$", ".log", s))
  t0 <- Sys.time()
  code <- system2(rscript, c("--encoding=UTF-8", shQuote(file.path("R", s))), stdout = log_file, stderr = log_file)
  status$seconds[i] <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  status$status[i] <- if (code == 0) "ok" else "error"
  message(sprintf("  %s (%.0f s)", if (code == 0) "done" else paste("FAILED - see", log_file), status$seconds[i]))
  if (code != 0 && stop_on_error) stop("Stopped at ", s)
}

write.csv(status, "output/logs/_run_status.csv", row.names = FALSE)
print(status, row.names = FALSE)
message("Figures: output/figures/   tables: output/tables/   logs: output/logs/")
