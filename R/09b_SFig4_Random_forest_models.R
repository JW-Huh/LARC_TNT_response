## =============================================================================
##  09b_SFig4_Random_forest_models.R
##  ---------------------------------------------------------------------------
##  PANELS      Supplementary Figure 4B-J (results; the plots are drawn in 09c)
##
##  PURPOSE     Exploratory random-forest classification of pCR versus non-pCR
##              from each omic layer, with a global feature screen, a grid search
##              over screening cut-offs, patient-grouped 5-fold cross-validation
##              and SHAP feature attributions.
##
##  TWO MODES   refit_rf <- FALSE  (default)  replay: load the submitted model
##                  results (input/random_forest/rf_results_submitted.rds),
##                  re-verify every AUC from the stored out-of-fold predictions
##                  and pass them on.  Fast, exact reproduction of the figure.
##              refit_rf <- TRUE   full computation (~ 1-2 h): screening, grid
##                  search (3,600 / 900 cut-off combinations per assay), final
##                  cross-validated forests with importance and SHAP, and the
##                  clinical model.  Grid search, selected features and first-
##                  fold predictions match the submitted run; SHAP draws (and
##                  therefore later folds) differ slightly under the current
##                  package versions, so a refit is a NEW fit (see REVIEW notes).
##
##  INPUT       output/cache/random_forest_inputs.RData      (from 09a)
##              input/random_forest/rf_results_submitted.rds (replay mode)
##
##  OUTPUT      output/cache/random_forest_results.rds
##              output/tables/SFig4B_random_forest_AUC.csv
##              output/tables/SFig4_random_forest_selected_features.csv
##
##  METHODS  (randomForest 4.7-1.2, rsample 1.3.0, fastshap 0.1.1, pROC 1.19)
##   * Global screen (all 42 samples, both visits): for every feature, the
##     pCR - non-pCR pairwise differences at baseline and after RT give (i) the
##     direction concordance between visits, (ii) a binomial-tail score per
##     visit (geometric mean = "screening score"; not a calibrated P value
##     because samples enter many pairs), (iii) consistency = mean proportion of
##     pairs in the majority direction, (iv) D = 2 x AUC - 1 (pairs where both
##     values are 0 are ignored), (v) prevalence per group, (vi) a pooled
##     Wilcoxon P.
##   * Grid search: prevalence (0.1-0.5; 0.3-0.7 for KOs), |D| (0.1-0.3),
##     consistency (0.7-0.9), screening score (0.01-0.2; 0.01 only for
##     metabolites), Wilcoxon P (0.01-0.2) and top-N (10-50).  For each
##     combination the passing features (sorted by screening score, top N) feed
##     a 1,000-tree random forest evaluated by stratified 5-fold CV grouped by
##     patient (rsample::group_vfold_cv, seed 123); the pooled out-of-fold AUC
##     is recorded and the best combination (ties -> smaller N) is kept.
##     Screening and tuning use the same folds, so the AUC is optimistic.
##   * Final model: same CV with the selected features; per-fold MDA / MDG
##     converted to percentile ranks and averaged; SHAP values by fastshap
##     (100 Monte-Carlo simulations) on the held-out fold, mean |SHAP| per feature.
##   * Clinical model: cT stage code, cN stage and CEA >= 5 ng/mL, one row per
##     patient, same CV scheme.
##
##  R PACKAGES  pROC, dplyr, tidyr (+ randomForest, rsample, fastshap for a refit)
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
})
report_versions(c("pROC", "randomForest", "rsample", "fastshap"))

refit_rf <- FALSE          # <-- TRUE = recompute everything (slow); FALSE = replay submitted results

load(cache_file("random_forest_inputs.RData"))               # inputs, clinical (script 09a)

## Internal two-level labels used by the models ("CR" = pCR).  Kept as in the
## original code so that the seeded fold assignment is reproduced.
to_cr <- function(response) factor(ifelse(response == "pCR", "CR", "nonCR"), levels = c("nonCR", "CR"))
auc_of <- function(labels, prob) as.numeric(auc(roc(labels, prob, levels = c("nonCR", "CR"), direction = "<", quiet = TRUE)))


if (!refit_rf) {
  ## ---- A. Replay mode: verify and reuse the submitted results ---------------------
  all_results <- readRDS(input_file("random_forest", "rf_results_submitted.rds"))
  stopifnot(all(c(names(inputs), "CI") %in% names(all_results)))
  for (nm in names(inputs)) {
    res <- all_results[[nm]]; x <- inputs[[nm]]
    stopifnot(
      identical(res$Pred$SampleID, rownames(x$values)),
      identical(as.character(res$Pred$TRG_1), as.character(to_cr(x$metadata$Response))),
      all(res$Selected_features %in% colnames(x$values))
    )
    stopifnot(abs(auc_of(res$Pred$TRG_1, res$Pred$Pred_CR) - res$AUC) < 1e-12)   # AUC recomputed from the predictions
  }
  stopifnot(abs(auc_of(all_results$CI$Pred$TRG_1, all_results$CI$Pred$Pred_CR) - all_results$CI$AUC) < 1e-12)
  message("Replay mode: submitted out-of-fold predictions verified; no model was refitted.")

} else {
  ## ---- B. Full refit ------------------------------------------------------------------
  suppressPackageStartupMessages({ library(randomForest); library(rsample); library(fastshap) })

  SEED <- 123; N_TREE <- 1000; V_FOLD <- 5; NSIM_SHAP <- 100      # <-- model settings
  shap_fun <- function(object, newdata) predict(object, newdata, type = "prob")[, "CR"]

  ## B1. Global screening statistics for one assay ---------------------------------
  screen_features <- function(values, meta) {
    cr  <- unique(meta$SubjectID[meta$Response == "pCR"])
    ncr <- unique(meta$SubjectID[meta$Response == "non-pCR"])
    b <- values[meta$Timepoint == "Before", , drop = FALSE];  rownames(b) <- meta$SubjectID[meta$Timepoint == "Before"]
    o <- values[meta$Timepoint == "Ongoing", , drop = FALSE]; rownames(o) <- meta$SubjectID[meta$Timepoint == "Ongoing"]
    cr_b <- intersect(cr, rownames(b)); ncr_b <- intersect(ncr, rownames(b))
    cr_o <- intersect(cr, rownames(o)); ncr_o <- intersect(ncr, rownames(o))

    one_feature <- function(f) {
      b_d <- as.vector(outer(b[cr_b, f], b[ncr_b, f], "-"))       # all pCR - non-pCR pairs, baseline
      o_d <- as.vector(outer(o[cr_o, f], o[ncr_o, f], "-"))       # ... after RT
      b_pos <- sum(b_d > 0); b_neg <- sum(b_d < 0); b_n <- b_pos + b_neg
      o_pos <- sum(o_d > 0); o_neg <- sum(o_d < 0); o_n <- o_pos + o_neg
      all_cr <- c(b[cr_b, f], o[cr_o, f]); all_ncr <- c(b[ncr_b, f], o[ncr_o, f])

      ## D = 2 * AUC - 1 over all pooled pairs (pairs with both values 0 ignored)
      gt <- outer(all_cr, all_ncr, ">"); eq <- outer(all_cr, all_ncr, "==")
      both0 <- outer(all_cr == 0, all_ncr == 0, "&"); valid <- !both0
      auc_f <- if (sum(valid) > 0) (sum(gt & valid) + 0.5 * sum(eq & valid)) / sum(valid) else NA_real_
      D <- if (is.na(auc_f)) NA_real_ else 2 * auc_f - 1

      wp <- tryCatch(suppressWarnings(wilcox.test(all_cr, all_ncr, exact = FALSE)$p.value), error = function(e) NA_real_)

      ## direction concordance between visits -> consistency and binomial-tail score
      cons <- pval <- NA_real_
      if (b_n > 0 && o_n > 0) {
        b_dir <- ifelse(b_pos >= b_neg, 1L, -1L); o_dir <- ifelse(o_pos >= o_neg, 1L, -1L)
        if (b_dir == o_dir) {
          b_pv <- binom.test(max(b_pos, b_neg), b_n, 0.5, "greater")$p.value
          o_pv <- binom.test(max(o_pos, o_neg), o_n, 0.5, "greater")$p.value
          pval <- sqrt(b_pv * o_pv)                              # geometric mean = screening score
          cons <- (max(b_pos, b_neg) / b_n + max(o_pos, o_neg) / o_n) / 2
        }
      }
      data.frame(Feature = f, Prev_CR = mean(all_cr > 0), Prev_nonCR = mean(all_ncr > 0),
                 D = D, Consistency = cons, Pvalue = pval, Wilcoxon_p = wp)
    }
    bind_rows(lapply(colnames(values), one_feature))
  }

  ## B2. Apply one set of cut-offs to the screening table ----------------------------
  select_features <- function(scores, g) {
    sel <- scores %>%
      filter(!is.na(D), !is.na(Consistency), !is.na(Pvalue)) %>%
      filter(Prev_CR >= g$prev_cutoff | Prev_nonCR >= g$prev_cutoff) %>%
      filter(abs(D) >= g$d_cutoff, Consistency >= g$cons_cutoff, Pvalue < g$pval_cutoff) %>%
      filter(!is.na(Wilcoxon_p), Wilcoxon_p < g$wilcox_cutoff) %>%
      arrange(Pvalue)
    if (nrow(sel) > g$n_top) sel <- slice_head(sel, n = g$n_top)
    sel
  }

  ## B3. Patient-grouped, stratified 5-fold CV of a random forest ----------------------
  ## Returns pooled out-of-fold P(CR) and, when `full = TRUE`, importance and SHAP.
  cv_forest <- function(rf_data, features, folds, full = FALSE) {
    pred <- setNames(rep(NA_real_, nrow(rf_data)), rownames(rf_data))
    imp <- shap <- vector("list", nrow(folds))
    set.seed(SEED)
    for (i in seq_len(nrow(folds))) {
      tr <- analysis(folds$splits[[i]]); te <- assessment(folds$splits[[i]])
      mod <- randomForest(x = tr[, features, drop = FALSE], y = tr$TRG_1, ntree = N_TREE, importance = full)
      pred[match(rownames(te), rownames(rf_data))] <- predict(mod, te[, features, drop = FALSE], type = "prob")[, "CR"]
      if (full) {
        imp[[i]] <- as.data.frame(importance(mod)) %>% tibble::rownames_to_column("Feature") %>%
          transmute(Feature, MDA = MeanDecreaseAccuracy, MDG = MeanDecreaseGini, Fold = i)
        sv <- explain(mod, X = tr[, features, drop = FALSE], newdata = te[, features, drop = FALSE],
                      pred_wrapper = shap_fun, nsim = NSIM_SHAP)
        shap[[i]] <- as.data.frame(sv) %>% mutate(SampleID = rownames(te)) %>%
          pivot_longer(-SampleID, names_to = "Feature", values_to = "SHAP") %>%
          left_join(as.data.frame(te[, features, drop = FALSE]) %>% mutate(SampleID = rownames(te)) %>%
                      pivot_longer(-SampleID, names_to = "Feature", values_to = "Abundance"),
                    by = c("SampleID", "Feature")) %>%
          mutate(Fold = i)
        cat(sprintf("    fold %d/%d done\n", i, nrow(folds)))
      }
    }
    list(pred = pred, importance = bind_rows(imp), shap = bind_rows(shap))
  }

  ## B4. One assay from screening to final model ----------------------------------------
  all_results <- list()
  for (nm in names(inputs)) {
    message("Random forest: ", nm)
    values <- inputs[[nm]]$values; meta <- inputs[[nm]]$metadata
    scores <- screen_features(values, meta)

    ## assay-specific cut-off grid (as recorded in the submitted results)
    grid <- expand.grid(
      prev_cutoff   = if (nm == "KO") c(0.3, 0.4, 0.5, 0.6, 0.7) else c(0.1, 0.2, 0.3, 0.4, 0.5),
      d_cutoff      = c(0.1, 0.2, 0.3),
      cons_cutoff   = c(0.7, 0.8, 0.9),
      pval_cutoff   = if (nm == "Metabolite") 0.01 else c(0.01, 0.05, 0.1, 0.2),
      wilcox_cutoff = c(0.01, 0.05, 0.1, 0.2),
      n_top         = c(10, 20, 30, 40, 50)
    )

    rf_data <- as.data.frame(values)
    rf_data$SNU_ID <- meta$SubjectID                     # grouping variable (all visits of a patient together)
    rf_data$TRG_1  <- to_cr(meta$Response)
    set.seed(SEED)
    folds <- group_vfold_cv(rf_data, group = SNU_ID, v = V_FOLD, strata = TRG_1)
    labels <- rf_data$TRG_1

    ## grid search (identical feature sets are evaluated once)
    grid$AUC <- NA_real_; seen <- new.env(parent = emptyenv())
    for (gi in seq_len(nrow(grid))) {
      fc <- select_features(scores, grid[gi, ])$Feature
      if (length(fc) == 0) next
      key <- paste(fc, collapse = "\034")
      if (exists(key, seen, inherits = FALSE)) { grid$AUC[gi] <- get(key, seen); next }
      pred <- cv_forest(rf_data, fc, folds)$pred
      ok <- !is.na(pred)
      if (length(unique(labels[ok])) == 2) grid$AUC[gi] <- tryCatch(auc_of(labels[ok], pred[ok]), error = function(e) NA_real_)
      assign(key, grid$AUC[gi], seen)
    }
    stopifnot(any(is.finite(grid$AUC)))
    best <- grid %>% filter(!is.na(AUC)) %>% arrange(desc(AUC), n_top) %>% slice(1)
    cat("  best cut-offs:\n"); print(best, row.names = FALSE)

    ## final model with importance and SHAP
    sel <- select_features(scores, best); fc <- sel$Feature
    cat(sprintf("  %d features selected\n", length(fc)))
    fit <- cv_forest(rf_data, fc, folds, full = TRUE)
    ok <- !is.na(fit$pred)
    roc_obj <- roc(labels[ok], fit$pred[ok], levels = c("nonCR", "CR"), direction = "<", quiet = TRUE)

    shap_summary <- fit$shap %>% group_by(Feature) %>% summarise(MeanAbsSHAP = mean(abs(SHAP)), .groups = "drop")
    importance_summary <- fit$importance %>% group_by(Fold) %>%
      mutate(MDA = percent_rank(MDA), MDG = percent_rank(MDG)) %>% ungroup() %>%
      group_by(Feature) %>% summarise(MDA = mean(MDA), MDG = mean(MDG), .groups = "drop")
    feature_summary <- shap_summary %>% left_join(importance_summary, by = "Feature") %>%
      left_join(select(sel, Feature, Consistency, Pvalue), by = "Feature") %>% arrange(desc(MeanAbsSHAP))

    all_results[[nm]] <- list(
      Input = nm, Model = "RF_Global", FS = "Direction", Grid = grid, Best_criteria = best,
      Scores_global = scores, Selected_features = fc, Selected_scores = sel,
      AUC = as.numeric(auc(roc_obj)), ROC = roc_obj,
      Pred = data.frame(SampleID = rownames(rf_data), SubjectID = rf_data$SNU_ID, TRG_1 = labels, Pred_CR = fit$pred),
      Importance = fit$importance, Importance_summary = importance_summary,
      SHAP = fit$shap, SHAP_summary = shap_summary, Feature_summary = feature_summary
    )
    cat(sprintf("  pooled out-of-fold AUC = %.3f\n", all_results[[nm]]$AUC))
  }

  ## B5. Clinical model (one row per patient) ------------------------------------------
  ci_features <- c("cT_stage_code", "cN_stage", "CEA_high")
  rf_ci <- clinical %>% transmute(SNU_ID = SubjectID, TRG_1 = to_cr(Response), across(all_of(ci_features)))
  rownames(rf_ci) <- rf_ci$SNU_ID
  set.seed(SEED)
  folds_ci <- group_vfold_cv(rf_ci, group = SNU_ID, v = V_FOLD, strata = TRG_1)
  fit_ci <- cv_forest(rf_ci, ci_features, folds_ci)
  roc_ci <- roc(rf_ci$TRG_1, fit_ci$pred, levels = c("nonCR", "CR"), direction = "<", quiet = TRUE)
  all_results[["CI"]] <- list(Input = "CI", Model = "RF_CI", FS = NA, AUC = as.numeric(auc(roc_ci)), ROC = roc_ci,
                              Pred = data.frame(SubjectID = rf_ci$SNU_ID, TRG_1 = rf_ci$TRG_1, Pred_CR = fit_ci$pred))
  cat(sprintf("Clinical model: pooled out-of-fold AUC = %.3f\n", all_results$CI$AUC))
}


## ---- C. Save results and summary tables ---------------------------------------------

saveRDS(all_results, cache_file("random_forest_results.rds"))

auc_tbl <- data.frame(Input = names(all_results),
                      AUC = vapply(all_results, `[[`, numeric(1), "AUC"),
                      N_samples = vapply(all_results, function(r) nrow(r$Pred), integer(1)),
                      N_features = vapply(all_results, function(r) if (is.null(r$Selected_features)) 3L else length(r$Selected_features), integer(1)),
                      Mode = ifelse(refit_rf, "refit", "submitted"), row.names = NULL)
print(auc_tbl, row.names = FALSE, digits = 4)
save_table(auc_tbl, "SFig4B_random_forest_AUC")

selected_tbl <- bind_rows(lapply(names(inputs), function(nm) {
  cbind(Input = nm, all_results[[nm]]$Feature_summary)
}))
save_table(selected_tbl, "SFig4_random_forest_selected_features")

message("Done: 09b_SFig4_Random_forest_models.R")
