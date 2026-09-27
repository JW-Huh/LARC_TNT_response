## =============================================================================
##  20_Fig5_MOFA_preparation_and_training.R
##  ---------------------------------------------------------------------------
##  PANELS      (no figure) - MOFA input matrices and the trained factor model
##              used by scripts 21-23 (Figure 5B-F, Supplementary Figure 6)
##
##  PURPOSE     Build the four MOFA views from the preprocessed blocks of script
##              18, check them against the frozen model input, and load the
##              frozen trained model (factor scores, weights, variance).
##
##  TWO MODES   refit_mofa <- FALSE (default)  rebuild the views, verify that they
##                  equal the frozen training input (relative tolerance 1e-4) and use the
##                  frozen trained model.  MOFA training is stochastic and
##                  depends on the Python back end, so the frozen model is the
##                  record of the submitted analysis.
##              refit_mofa <- TRUE   train the 4 x 2 model grid (K = 6/10/14 with
##                  spike-slab weights, K = 10 without; seeds 41 and 97) with
##                  MOFA2 and select the model with the highest ELBO.  Requires
##                  MOFA2 and a mofapy2 Python environment reachable by
##                  reticulate (set RETICULATE_PYTHON before starting R).  A
##                  refit can change factor numbering and sign.
##
##  INPUT       output/cache/four_omics_inputs.RData          (script 18)
##              input/mofa/mofa_input_frozen.rds              frozen training input
##              input/mofa/mofa_trained_model_frozen.rds      frozen trained model tables
##              input/metadata/patients.csv
##
##  OUTPUT      output/cache/mofa_model.RData   (mofa_input, mofa_model)
##              output/tables/Fig5B_MOFA_training_log.csv, Fig5B_MOFA_input_dimensions.csv
##
##  METHODS     MOFA2 1.12.1 / mofapy2: Bayesian group factor analysis with
##              Gaussian likelihoods, ARD on weights, spike-slab weights, view
##              scaling, 42 samples (all metagenome samples; absent views are
##              missing data), up to 2,000 features per view.  Views: species
##              (Hellinger, 392), KO (Hellinger, 2,000), metabolites (floor-
##              corrected, z-scored, 25), host (VST, 2,000).  The selected model
##              has K = 14 factors (spike-slab, seed 41).
##
##  R PACKAGES  (MOFA2, reticulate for a refit)
## =============================================================================

source("R/_common.R")
source("R/_helpers_mofa.R")
report_versions(c("MOFA2", "nlme"))

refit_mofa <- FALSE          # <-- TRUE = retrain (see header)
view_tolerance <- 1e-4       # relative tolerance for "rebuilt view == frozen training input"
                             # (species, metabolite and host views are identical; the KO view
                             #  differs by 1.5e-5 relative / < 5e-6 absolute, i.e. a marginal
                             #  KO in the sample totals - no feature or ranking changes)

load(cache_file("four_omics_inputs.RData"))                          # coherence (script 18)
patients <- read_patients()


## ---- 1. Rebuild the MOFA views ---------------------------------------------------------

mofa_input <- build_mofa_input(coherence, patients)
dims <- data.frame(view = names(mofa_input$values), samples = sapply(mofa_input$values, nrow),
                   features = sapply(mofa_input$values, ncol),
                   observed_samples = sapply(mofa_input$values, function(a) sum(rowSums(is.finite(a)) > 0)), row.names = NULL)
print(dims, row.names = FALSE)
save_table(dims, "Fig5B_MOFA_input_dimensions")


## ---- 2. Frozen model (default) or refit --------------------------------------------------

if (!refit_mofa) {
  frozen <- readRDS(input_file("mofa", "mofa_input_frozen.rds"))
  ## The laboratory table used for the submission spelled one metabolite
  ## "Lithocholi_acid"; the curated input uses the correct name, so the frozen
  ## objects are renamed on loading (values only - nothing else changes).
  rename_frozen <- function(x) { x[x == "Lithocholi_acid"] <- "Lithocholic_acid"; x }
  colnames(frozen$values$metabolite) <- rename_frozen(colnames(frozen$values$metabolite))
  mismatch <- character()
  for (v in names(mofa_input$values)) {
    a <- mofa_input$values[[v]]; b <- frozen$values[[v]]
    rr <- intersect(rownames(a), rownames(b)); cc <- intersect(colnames(a), colnames(b))
    d <- suppressWarnings(max(abs(a[rr, cc, drop = FALSE] - b[rr, cc, drop = FALSE]), na.rm = TRUE))
    cat(sprintf("%-10s rebuilt %d x %d | frozen %d x %d | shared %d x %d | max |diff| %.2e | same NA pattern %s\n",
                v, nrow(a), ncol(a), nrow(b), ncol(b), length(rr), length(cc), d,
                identical(is.na(a[rr, cc, drop = FALSE]), is.na(b[rr, cc, drop = FALSE]))))
    ## features are compared by name: the order of the 2,000 top-variance
    ## features can differ between runs (ties), which is irrelevant for MOFA
    if (setequal(rownames(a), rownames(b)) && setequal(colnames(a), colnames(b))) b <- b[rownames(a), colnames(a), drop = FALSE]
    eq <- all.equal(a, b, tolerance = view_tolerance, check.attributes = TRUE)
    if (!isTRUE(eq)) {
      mismatch <- c(mismatch, v); cat("   ", paste(eq, collapse = "; "), "\n")
      cat("   only rebuilt:", paste(head(setdiff(colnames(a), colnames(b)), 5), collapse = ", "),
          "| only frozen:", paste(head(setdiff(colnames(b), colnames(a)), 5), collapse = ", "), "\n")
    }
  }
  if (length(mismatch)) stop("Rebuilt MOFA view(s) ", paste(mismatch, collapse = ", "),
                             " differ from the frozen training input; inspect before reusing the frozen factors.")
  stopifnot(identical(mofa_input$metadata$SampleID, frozen$metadata$SampleID))
  message("All four rebuilt views match the frozen training input (relative tolerance ", view_tolerance, ").")
  mofa_model <- readRDS(input_file("mofa", "mofa_trained_model_frozen.rds"))
  mofa_model$extracted$weights$feature <- rename_frozen(mofa_model$extracted$weights$feature)
  if ("feature_label" %in% names(mofa_model$extracted$weights))
    mofa_model$extracted$weights$feature_label <- rename_frozen(mofa_model$extracted$weights$feature_label)
} else {
  mofa_model <- train_mofa_grid(mofa_input$values, mofa_input$metadata, file.path(paths$cache, "mofa_models"))
  message("Refit finished. Factor numbering / sign may differ from the submission - check `display_factors` in R/_helpers_mofa.R.")
}

print(mofa_model$training_log[, c("K", "spikeslab", "seed", "ELBO", "selected")], row.names = FALSE)
save_table(mofa_model$training_log[, c("K", "spikeslab", "seed", "ELBO", "selected")], "Fig5B_MOFA_training_log")
save(mofa_input, mofa_model, file = cache_file("mofa_model.RData"))

message("Done: 20_Fig5_MOFA_preparation_and_training.R")
