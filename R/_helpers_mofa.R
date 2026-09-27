## =============================================================================
##  R/_helpers_mofa.R  --  functions shared by the MOFA scripts 20-23
##  ---------------------------------------------------------------------------
##  Sourced after R/_common.R.  Nothing here runs an analysis by itself.
##
##  Contents
##    mofa_cfg                 all tuning constants of the MOFA analysis
##    build_mofa_input()       four preprocessed views + sample metadata (script 20)
##    train_mofa_grid()        optional refit with MOFA2 (script 20, refit_mofa = TRUE)
##    hedges(), patient_means(), permutation_labels(), fit_mixed(), gls_stat() ...
##                             statistics used by analyze_association()
##    analyze_association()    factor scores ~ response (script 21)
##    period_effects()         Hedges' g overall / baseline / after RT (script 21)
##    variance_explained()     reconstruct per-view and total R2 from scores x weights
##    candidate_features()     which features to display for a factor (script 22)
##    partial_spearman()       covariate-adjusted rank correlation with permutations
##    host_rank_screen()       Pearson / Spearman / partial Spearman for TJP1 (script 23)
## =============================================================================

suppressPackageStartupMessages(library(nlme))

## ---- settings ----------------------------------------------------------------------

mofa_cfg <- list(
  seed = 20260908L, permutations = 999L, bootstrap = 999L,
  minimum_observed_views = 1L, minimum_samples = 16L,
  prevalence = 0.20,                          # species / KO prevalence filter
  minimum_feature_observations = 4L,
  feature_caps = c(species = 2000L, ko = 2000L, metabolite = 2000L, host = 2000L),   # top-variance features per view
  metabolite_scale_features = TRUE,
  metabolite_floor_min_repeats = 2L, metabolite_floor_relative_tolerance = 1e-8, metabolite_floor_absolute_tolerance = 1e-12,
  metabolite_finite_fraction = 0.50, metabolite_detected_fraction = 0.50,
  display_factors = c("Factor7", "Factor11"),
  stratify_permutations_by_visit_pattern = TRUE,
  minimum_view_r2 = 0.01, minimum_normalized_loading = 0.25, minimum_display_loading = 0.10,
  minimum_correlation_patients = 8L,
  training_grid = data.frame(K = c(6L, 10L, 14L, 10L), spikeslab = c(TRUE, TRUE, TRUE, FALSE)),
  training_seeds = c(41L, 97L), maxiter = 2000L
)
view_order <- c("species", "ko", "metabolite", "host")


## ---- MOFA input ---------------------------------------------------------------------

## Preprocess the four views for MOFA from the blocks of script 18.
##   species / KO : prevalence >= 20 % (of the 42 metagenomes), relative
##                  abundance renormalised over the retained features (the
##                  block's `relative` matrix of script 18) -> Hellinger (sqrt)
##   metabolite   : repeated detection-floor values set to 0; >= 50 % finite and
##                  >= 50 % detected; centred and scaled
##   host         : blind DESeq2 VST (all count-QC genes)
##   all views    : non-constant features observed in >= 4 samples, top 2,000 by
##                  variance, centred (features x samples arranged on the 42
##                  metagenome samples; unobserved views are NA)
## NOTE  the sample totals are taken AFTER the prevalence filter (as in the
##       submitted run); dividing by the total of all features instead changes
##       the Hellinger values by up to 0.1 and the frozen model no longer matches.
build_mofa_input <- function(coherence, patients) {
  cfg <- mofa_cfg
  mods <- coherence$modalities
  ## sample universe: every SubjectID__Timepoint observed in any view
  meta <- unique(do.call(rbind, lapply(mods, function(b) b$sample_meta[, c("SampleID", "SubjectID", "Timepoint")])))
  meta <- meta[order(meta$SubjectID, meta$Timepoint), ]
  values <- list()
  for (v in view_order) {
    b <- mods[[v]]
    if (v %in% c("species", "ko")) {
      a <- b$raw; a <- a[rowSums(a) > 0, , drop = FALSE]
      keep <- colMeans(a > 0) >= cfg$prevalence & colSums(a) > 0
      a <- a[, keep, drop = FALSE]
      a <- sqrt(a / rowSums(a))                                                    # totals over the retained features
    } else if (v == "metabolite") {
      a <- b$raw
      for (j in seq_len(ncol(a))) {
        x <- a[, j]; floor <- min(x, na.rm = TRUE)
        tol <- max(cfg$metabolite_floor_absolute_tolerance, abs(floor) * cfg$metabolite_floor_relative_tolerance)
        at_floor <- is.finite(x) & abs(x - floor) <= tol
        if (sum(at_floor) >= cfg$metabolite_floor_min_repeats) a[at_floor, j] <- 0     # repeated floor -> 0
      }
      a <- a[rowMeans(is.finite(a)) >= cfg$metabolite_finite_fraction, , drop = FALSE]
      keep <- colMeans(is.finite(a)) >= cfg$metabolite_finite_fraction & colMeans(is.finite(a) & a > 0) >= cfg$metabolite_detected_fraction
      a <- a[, keep, drop = FALSE]
    } else {
      a <- b$vst_qc                                                                 # biopsies x all QC genes
    }
    vv <- apply(a, 2, var, na.rm = TRUE)
    keep <- which(is.finite(vv) & vv > 1e-12 & colSums(is.finite(a)) >= cfg$minimum_feature_observations)
    keep <- head(keep[order(-vv[keep], colnames(a)[keep])], cfg$feature_caps[[v]])
    a <- a[, keep, drop = FALSE]
    a <- sweep(a, 2, colMeans(a, na.rm = TRUE), "-")
    if (v == "metabolite" && cfg$metabolite_scale_features) a <- sweep(a, 2, apply(a, 2, sd, na.rm = TRUE), "/")
    m <- matrix(NA_real_, nrow(meta), ncol(a), dimnames = list(meta$SampleID, colnames(a)))
    m[match(rownames(a), meta$SampleID), ] <- a
    values[[v]] <- m
  }
  observed <- sapply(values, function(a) rowSums(is.finite(a)) > 0)
  keep <- rowSums(observed) >= cfg$minimum_observed_views
  meta <- meta[keep, ]; values <- lapply(values, function(a) a[keep, , drop = FALSE])
  meta$TRG_plot <- ifelse(patients$Response[match(meta$SubjectID, patients$SubjectID)] == "pCR", "pCR", "non_pCR")
  meta$TRG_score <- patients$TRG_score[match(meta$SubjectID, patients$SubjectID)]
  meta$Response <- as.integer(meta$TRG_plot == "pCR")
  meta$Time <- as.integer(meta$Timepoint == "Ongoing")
  meta$n_observed_views <- rowSums(observed[keep, ])
  rownames(meta) <- NULL
  list(values = values, metadata = meta)
}

## Optional refit: grid of K x spike-slab settings x seeds, best ELBO wins.
## Requires MOFA2 and a working mofapy2 Python environment (see script 20).
train_mofa_grid <- function(values, meta, model_dir) {
  cfg <- mofa_cfg
  suppressPackageStartupMessages(library(MOFA2))
  input <- lapply(values, function(a) {
    keep <- colSums(is.finite(a)) >= 4 & apply(a, 2, var, na.rm = TRUE) > 1e-12
    t(a[, keep, drop = FALSE])
  })
  grid <- unique(transform(cfg$training_grid, K = pmin(K, nrow(meta) - 2L)))
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  log <- data.frame()
  for (i in seq_len(nrow(grid))) for (seed in cfg$training_seeds) {
    tag <- sprintf("pooled_K%d_spike%d_seed%d", grid$K[i], as.integer(grid$spikeslab[i]), seed)
    path <- file.path(model_dir, paste0(tag, ".hdf5"))
    message("Training ", tag)
    if (!file.exists(path)) {
      obj <- create_mofa(input)
      dop <- get_default_data_options(obj); dop$scale_views <- TRUE; dop$scale_groups <- FALSE; dop$center_groups <- TRUE
      mop <- get_default_model_options(obj); mop$num_factors <- grid$K[i]; mop$likelihoods[] <- "gaussian"
      mop$ard_weights <- TRUE; mop$ard_factors <- FALSE; mop$spikeslab_weights <- grid$spikeslab[i]; mop$spikeslab_factors <- FALSE
      top <- get_default_training_options(obj); top$maxiter <- cfg$maxiter; top$convergence_mode <- "slow"; top$seed <- seed
      top$drop_factor_threshold <- -1; top$startELBO <- 1L; top$freqELBO <- 5L; top$verbose <- FALSE
      obj <- prepare_mofa(obj, data_options = dop, model_options = mop, training_options = top)
      run_mofa(obj, outfile = path, save_data = TRUE, use_basilisk = FALSE)
    }
    fit <- load_model(path)
    log <- rbind(log, data.frame(K = grid$K[i], spikeslab = grid$spikeslab[i], seed = seed, ELBO = as.numeric(get_elbo(fit)), path = path))
  }
  log$selected <- seq_len(nrow(log)) == which.max(log$ELBO)
  model <- load_model(log$path[log$selected])
  z <- as.matrix(get_factors(model)[[1]])
  scores <- data.frame(meta, z[match(meta$SampleID, rownames(z)), ], check.names = FALSE)
  ww <- get_weights(model)
  weights <- do.call(rbind, lapply(names(ww), function(v) data.frame(view = v, feature = rep(rownames(ww[[v]]), ncol(ww[[v]])),
                                                                     factor = rep(colnames(ww[[v]]), each = nrow(ww[[v]])), weight = as.vector(ww[[v]]))))
  weights$feature_label <- weights$feature
  variance <- as.data.frame(get_variance_explained(model, as.data.frame = TRUE)$r2_per_factor)
  variance$r2 <- variance$value / 100
  list(extracted = list(scores = scores, weights = weights, variance = variance), training_log = log, selected_hdf5 = log$path[log$selected])
}


## ---- small statistics ------------------------------------------------------------------

## Hedges' g of x between y == 1 (pCR) and y == 0 (non-pCR)
hedges <- function(x, y) {
  a <- x[is.finite(y) & y == 1 & is.finite(x)]; b <- x[is.finite(y) & y == 0 & is.finite(x)]
  if (min(length(a), length(b)) < 2) return(NA_real_)
  s <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  if (!is.finite(s) || s <= 0) return(NA_real_)
  (1 - 3 / (4 * (length(a) + length(b)) - 9)) * (mean(a) - mean(b)) / s
}

## Mean of each column per patient (rows = samples of d$SubjectID)
patient_means <- function(x, d) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  ids <- unique(d$SubjectID)
  out <- t(sapply(ids, function(id) colMeans(x[d$SubjectID == id, , drop = FALSE], na.rm = TRUE)))
  if (ncol(x) == 1) out <- matrix(out, ncol = 1)
  out[!is.finite(out)] <- NA; dimnames(out) <- list(ids, colnames(x)); out
}

## B permutations of the patient-level response label; when `restricted`,
## labels are only exchanged among patients with the same visit pattern
permutation_labels <- function(d, B, seed, restricted = TRUE) {
  ids <- unique(d$SubjectID); y <- setNames(d$Response[match(ids, d$SubjectID)], ids)
  pattern <- vapply(ids, function(id) paste(sort(d$Time[d$SubjectID == id]), collapse = ":"), "")
  blocks <- if (restricted) split(seq_along(ids), pattern) else list(seq_along(ids))
  set.seed(seed)
  out <- matrix(rep(y, each = B), B, length(ids), dimnames = list(NULL, ids))
  for (b in seq_len(B)) for (jj in blocks) out[b, jj] <- y[jj[sample.int(length(jj))]]
  out
}

## Random-intercept model (nlme) of a factor score on response
fit_mixed <- function(value, d, kind = c("pooled", "null")) {
  kind <- match.arg(kind)
  a <- data.frame(Value = value, Response = d$Response, SubjectID = d$SubjectID); a <- a[is.finite(a$Value), ]
  if (nrow(a) < 8 || length(unique(a$SubjectID)) < 6) return(NULL)
  form <- if (kind == "pooled") Value ~ Response else Value ~ 1
  for (opt in c("nlminb", "optim")) {
    f <- tryCatch(if (anyDuplicated(a$SubjectID)) lme(form, random = ~ 1 | SubjectID, data = a, method = "REML",
                                                        control = lmeControl(opt = opt, msMaxIter = 200, returnObject = FALSE))
                  else gls(form, data = a, method = "REML"), error = function(e) NULL)
    if (!is.null(f)) return(f)
  }
  NULL
}
coefficient <- function(f, term = "Response") {
  na <- c(beta = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_, P = NA_real_)
  if (is.null(f)) return(na)
  tt <- summary(f)$tTable; if (!term %in% rownames(tt)) return(na)
  df <- if ("DF" %in% colnames(tt)) tt[term, "DF"] else nobs(f) - length(coef(f)); z <- qt(0.975, df)
  c(beta = tt[term, "Value"], se = tt[term, "Std.Error"], lower = tt[term, "Value"] - z * tt[term, "Std.Error"],
    upper = tt[term, "Value"] + z * tt[term, "Std.Error"], P = tt[term, "p-value"])
}

## Response-blind whitening matrix from the null random-intercept model, and
## the generalised-least-squares F statistic for response after whitening
null_whitener <- function(y, d) {
  fit <- fit_mixed(y, d, "null")
  if (!is.null(fit) && inherits(fit, "lme")) {
    vc <- VarCorr(fit)
    V <- as.numeric(vc[1, "Variance"]) * outer(d$SubjectID, d$SubjectID, "==") + diag(as.numeric(vc[nrow(vc), "Variance"]), nrow(d))
    return(list(L = t(chol(V)), status = "random-intercept covariance"))
  }
  list(L = diag(nrow(d)), status = "identity covariance")
}
gls_stat <- function(y, L, r) {
  yy <- forwardsolve(L, y); X <- forwardsolve(L, cbind(1, r)); Q <- qr(X)
  if (Q$rank < 2) return(c(F = 0, beta = NA_real_))
  ss <- sum(qr.resid(Q, yy)^2); ss0 <- sum(qr.resid(qr(forwardsolve(L, matrix(1, length(y), 1))), yy)^2)
  c(F = max(0, ss0 - ss) / max(ss / (length(y) - 2), .Machine$double.eps), beta = qr.coef(Q, yy)[2])
}

## PERMANOVA-type F and R2 of the response in a two-factor Euclidean plane
euclidean_pair <- function(Z, r, pairs) {
  if (length(unique(r)) < 2) return(list(F = rep(0, ncol(pairs)), R2 = rep(0, ncol(pairs))))
  total <- colSums(sweep(Z, 2, colMeans(Z))^2)
  within <- colSums(sweep(Z[r == 0, , drop = FALSE], 2, colMeans(Z[r == 0, , drop = FALSE]))^2) +
            colSums(sweep(Z[r == 1, , drop = FALSE], 2, colMeans(Z[r == 1, , drop = FALSE]))^2)
  bt <- pmax(0, total - within)
  B <- bt[pairs[1, ]] + bt[pairs[2, ]]; W <- within[pairs[1, ]] + within[pairs[2, ]]
  list(F = B / pmax(W, .Machine$double.eps) * (nrow(Z) - 2), R2 = B / pmax(B + W, .Machine$double.eps))
}
perm_p <- function(obs, null) (1 + colSums(sweep(null, 2, obs - 1e-12, ">="))) / (nrow(null) + 1)
adjust_max <- function(obs, null) vapply(obs, function(x) (1 + sum(apply(null, 1, max) >= x - 1e-12)) / (nrow(null) + 1), numeric(1))
forest_tier <- function(p) ifelse(!is.finite(p), "Not testable", ifelse(p < 0.05, "P < 0.05", ifelse(p < 0.10, "0.05 <= P < 0.10", "NS")))

## Sign that makes the displayed factor increase towards pCR
factor_sign <- function(a, k) {
  b <- a$association$beta[match(k, a$association$factor)]
  if (!is.finite(b)) b <- a$association$g[match(k, a$association$factor)]
  if (is.finite(b) && b < 0) -1 else 1
}


## ---- factor scores versus response ----------------------------------------------------

## For every factor: mixed-model beta, Hedges' g on patient means, a
## patient-label permutation P (GLS F statistic, response-blind covariance),
## BH q and max-F adjusted P, bootstrap CI of g; plus a PERMANOVA-type test for
## every factor pair (used for the Factor 7 x Factor 11 plane).
analyze_association <- function(extracted, d) {
  cfg <- mofa_cfg
  f <- grep("^Factor[0-9]+$", names(extracted$scores), value = TRUE)
  raw_all <- as.matrix(extracted$scores[match(d$SampleID, extracted$scores$SampleID), f]); rownames(raw_all) <- d$SampleID
  f <- f[apply(raw_all, 2, sd) > 1e-10]; raw_all <- raw_all[, f, drop = FALSE]
  Z_all <- scale(raw_all)
  use <- is.finite(d$Response); d <- d[use, ]; Z <- Z_all[use, , drop = FALSE]
  ids <- unique(d$SubjectID); y <- d$Response[match(ids, d$SubjectID)]
  U <- patient_means(Z, d)
  B <- cfg$permutations
  labels <- permutation_labels(d, B, cfg$seed, cfg$stratify_permutations_by_visit_pattern)
  row_index <- match(d$SubjectID, ids)
  whiteners <- lapply(f, function(k) null_whitener(Z[, k], d))
  F_obs <- vapply(seq_along(f), function(j) gls_stat(Z[, j], whiteners[[j]]$L, d$Response)[["F"]], numeric(1))
  F_null <- t(sapply(seq_len(B), function(b) vapply(seq_along(f), function(j) gls_stat(Z[, j], whiteners[[j]]$L, labels[b, row_index])[["F"]], numeric(1))))
  colnames(F_null) <- f

  assoc <- do.call(rbind, lapply(seq_along(f), function(j) {
    cf <- coefficient(fit_mixed(Z[, j], d, "pooled"))
    data.frame(factor = f[j], beta = cf[["beta"]], lower = cf[["lower"]], upper = cf[["upper"]], model_P = cf[["P"]],
               g = hedges(U[, j], y), Wilcoxon_P = suppressWarnings(wilcox.test(U[y == 1, j], U[y == 0, j], exact = FALSE)$p.value),
               n_pCR = sum(y == 1), n_non_pCR = sum(y == 0), n_observations = nrow(d), covariance = whiteners[[j]]$status)
  }))
  assoc$permutation_P <- perm_p(F_obs, F_null); assoc$permutation_q <- p.adjust(assoc$permutation_P, "BH"); assoc$maxF_P <- adjust_max(F_obs, F_null)

  set.seed(cfg$seed + 1L)
  boot <- t(replicate(B, {
    ii <- unlist(lapply(0:1, function(g) { a <- which(y == g); a[sample.int(length(a), replace = TRUE)] }))
    vapply(seq_along(f), function(j) hedges(U[ii, j], y[ii]), numeric(1))
  }))
  assoc$g_lower <- apply(boot, 2, quantile, 0.025, na.rm = TRUE); assoc$g_upper <- apply(boot, 2, quantile, 0.975, na.rm = TRUE)
  assoc$tier <- forest_tier(assoc$permutation_P)
  ord <- order(assoc$permutation_P, -abs(assoc$g), assoc$factor)
  axes <- head(unique(c(intersect(cfg$display_factors, f), assoc$factor[ord])), 2)

  pairs <- combn(seq_along(f), 2)
  po <- euclidean_pair(Z, d$Response, pairs)
  pn <- t(sapply(seq_len(B), function(b) euclidean_pair(Z, labels[b, row_index], pairs)$F))
  pair_results <- data.frame(factor_x = f[pairs[1, ]], factor_y = f[pairs[2, ]], F = po$F, R2 = po$R2, P = perm_p(po$F, pn), maxF_P = adjust_max(po$F, pn))
  best <- which((pair_results$factor_x == axes[1] & pair_results$factor_y == axes[2]) | (pair_results$factor_x == axes[2] & pair_results$factor_y == axes[1]))

  list(association = assoc, pair_results = pair_results, axes = axes, best_pair = best, metadata = d, Z = Z, Z_all = Z_all, U = U)
}

## Hedges' g per factor for all visits pooled (patient means), baseline only
## and after-RT only, with bootstrap CI and permutation P
period_effects <- function(a) {
  cfg <- mofa_cfg
  rows <- list()
  for (period in c("Overall", "Baseline", "After RT")) {
    ii <- switch(period, Overall = seq_len(nrow(a$metadata)), Baseline = which(a$metadata$Timepoint == "Before"), `After RT` = which(a$metadata$Timepoint == "Ongoing"))
    d <- a$metadata[ii, ]; z <- a$Z[ii, , drop = FALSE]
    u <- patient_means(z, d); y <- d$Response[match(rownames(u), d$SubjectID)]
    enough <- min(table(factor(y, levels = 0:1))) >= 3
    set.seed(cfg$seed + 620L)
    bi <- if (enough) replicate(cfg$bootstrap, unlist(lapply(0:1, function(g) { jj <- which(y == g); sample(jj, length(jj), replace = TRUE) }))) else NULL
    perm <- if (period != "Overall" && enough) permutation_labels(d, cfg$permutations, cfg$seed + 621L) else NULL
    for (k in colnames(z)) {
      boot <- if (enough) apply(bi, 2, function(j) hedges(u[j, k], y[j])) else NA
      ci <- if (sum(is.finite(boot)) >= 0.8 * cfg$bootstrap) quantile(boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE) else c(NA, NA)
      pv <- NA_real_
      if (period == "Overall") pv <- a$association$permutation_P[a$association$factor == k]
      else if (enough) {
        wh <- null_whitener(z[, k], d); ob <- gls_stat(z[, k], wh$L, d$Response)[["F"]]
        np <- apply(perm, 1, function(r) gls_stat(z[, k], wh$L, r[match(d$SubjectID, colnames(perm))])[["F"]])
        pv <- (1 + sum(np >= ob - 1e-12)) / (length(np) + 1)
      }
      rows[[length(rows) + 1]] <- data.frame(factor = k, period = period, g = hedges(u[, k], y), lower = ci[1], upper = ci[2], P = pv,
                                             n_pCR = sum(y == 1), n_non_pCR = sum(y == 0))
    }
  }
  do.call(rbind, rows)
}

## Per-view and total variance explained, reconstructed from scores x weights
## on the model-scaled data (each view centred and divided by its overall SD)
variance_explained <- function(extracted, values) {
  factors <- unique(extracted$weights$factor)
  comp <- list()
  for (v in names(values)) {
    ww <- extracted$weights[extracted$weights$view == v, ]
    features <- unique(ww$feature)
    Y <- as.matrix(values[[v]])[, features, drop = FALSE]
    Z <- as.matrix(extracted$scores[match(rownames(Y), extracted$scores$SampleID), factors])
    Y <- sweep(Y, 2, colMeans(Y, na.rm = TRUE)); Y <- Y / sqrt(mean((Y - mean(Y, na.rm = TRUE))^2, na.rm = TRUE))
    W <- matrix(NA_real_, length(features), length(factors), dimnames = list(features, factors))
    W[cbind(match(ww$feature, features), match(ww$factor, factors))] <- ww$weight
    SST <- sum(Y^2, na.rm = TRUE)
    for (k in factors) comp[[length(comp) + 1]] <- data.frame(view = v, factor = k, SST = SST, SSE = sum((Y - tcrossprod(Z[, k], W[, k]))^2, na.rm = TRUE))
    comp[[length(comp) + 1]] <- data.frame(view = v, factor = "ALL_FACTORS", SST = SST, SSE = sum((Y - tcrossprod(Z, W))^2, na.rm = TRUE))
  }
  comp <- do.call(rbind, comp); comp$r2 <- pmax(0, 1 - comp$SSE / comp$SST)
  total <- do.call(rbind, lapply(unique(comp$factor), function(k) {
    z <- comp[comp$factor == k, ]; data.frame(factor = k, r2 = max(0, 1 - sum(z$SSE) / sum(z$SST)))
  }))
  list(components = comp, total = total)
}


## ---- feature selection for one factor -----------------------------------------------

display_feature_label <- function(view, feature, label) {
  label <- as.character(label); bad <- is.na(label) | !nzchar(trimws(label)); label[bad] <- feature[bad]
  label <- gsub("_", " ", label, fixed = TRUE)
  label[view == "ko" & feature == "K06926"] <- "Putative ATPase-family protein"
  label <- gsub("Indolepropionic acid", "Indole propionic acid", label, ignore.case = TRUE)
  ko <- which(view == "ko"); for (i in ko) if (!grepl(feature[i], label[i], fixed = TRUE)) label[i] <- paste0(feature[i], ": ", label[i])
  sp <- which(view == "species"); label[sp] <- gsub("\\bsp[.]? +", "sp. ", label[sp], perl = TRUE); label[sp] <- gsub("AM49[ _-]+4BH", "AM49-4BH", label[sp])
  label
}

## Features shown for a factor: view R2 >= 1 %, within-view relative |weight|
## >= 0.25, factor-wide normalised |weight| > 0.10, at most five per view
candidate_features <- function(weights, variance, factor) {
  cfg <- mofa_cfg
  w <- weights[weights$factor == factor & is.finite(weights$weight), ]
  w$key <- paste(w$view, w$feature, sep = "::")
  if (!"feature_label" %in% names(w)) w$feature_label <- w$feature          # a refit has no annotation column
  w$feature_label <- display_feature_label(w$view, w$feature, w$feature_label)
  w$relative_weight <- ave(w$weight, w$view, FUN = function(x) if (max(abs(x)) > 0) x / max(abs(x)) else 0 * x)
  w$factor_scaled_weight <- w$weight / max(abs(w$weight))
  w$r2 <- variance$r2[match(paste(w$view, w$factor), paste(variance$view, variance$factor))]
  w <- w[order(-abs(w$weight), w$key), ]
  active <- w[is.finite(w$r2) & w$r2 >= cfg$minimum_view_r2 & abs(w$relative_weight) >= cfg$minimum_normalized_loading & abs(w$weight) > 1e-12, ]
  do.call(rbind, lapply(view_order, function(v) head(active[active$view == v & abs(active$factor_scaled_weight) > cfg$minimum_display_loading, ], 5)))
}


## ---- covariate-adjusted correlations ------------------------------------------------------

## Partial Spearman correlation of x and y given Age, BMI and Sex (rank residuals),
## with a residual-permutation P value (mofa_cfg$permutations)
partial_spearman <- function(x, y, covariates, seed, adjust = c("Age", "BMI", "Sex")) {
  cfg <- mofa_cfg
  out <- list(n = 0L, rho = NA_real_, P = NA_real_, t_P = NA_real_, df = NA_integer_)
  ok <- is.finite(x) & is.finite(y) & complete.cases(covariates[, adjust, drop = FALSE])
  x <- x[ok]; y <- y[ok]; cc <- covariates[ok, adjust, drop = FALSE]; n <- length(x); out$n <- n
  if (n < cfg$minimum_correlation_patients) return(out)
  C <- matrix(1, n, 1)
  for (nm in adjust) if (length(unique(cc[[nm]])) > 1) C <- cbind(C, if (nm == "Sex") model.matrix(~ factor(cc[[nm]]))[, -1, drop = FALSE] else rank(cc[[nm]]))
  Q <- qr(C); df <- n - Q$rank - 1L; out$df <- df
  if (df < 3) return(out)
  rx <- qr.resid(Q, rank(x)); ry <- qr.resid(Q, rank(y))
  if (sd(rx) < 1e-10 || sd(ry) < 1e-10) return(out)
  r <- cor(rx, ry); out$rho <- r; out$t_P <- 2 * pt(-abs(r) * sqrt(df / max(1 - r^2, 1e-14)), df)
  set.seed(seed)
  idx <- replicate(cfg$permutations, sample.int(n))
  rp <- qr.resid(Q, (rank(y) - ry) + matrix(ry[idx], nrow = n))
  nul <- abs(drop(crossprod(rx, rp))) / sqrt(sum(rx^2) * colSums(rp^2)); nul[!is.finite(nul)] <- 0
  out$P <- (1 + sum(nul >= abs(r) - 1e-12)) / (cfg$permutations + 1)
  out
}

## Pearson / Spearman / partial Spearman of a factor score with host genes
## (permutation P with 999 patient permutations of the score residuals)
host_rank_screen <- function(x, Y, cv, adjust = character(), method = c("spearman", "pearson")) {
  cfg <- mofa_cfg
  method <- match.arg(method); tf <- if (method == "spearman") rank else as.numeric
  Y <- as.matrix(Y)
  ok <- is.finite(x) & complete.cases(Y) & (if (length(adjust)) complete.cases(cv[, adjust, drop = FALSE]) else TRUE)
  C <- matrix(1, sum(ok), 1)
  for (nm in adjust) { v <- cv[[nm]][ok]; if (length(unique(v)) > 1) C <- cbind(C, if (nm == "Sex") model.matrix(~ factor(v))[, -1, drop = FALSE] else tf(v)) }
  Q <- qr(C); rx <- qr.resid(Q, tf(x[ok])); df <- sum(ok) - Q$rank - 1L
  set.seed(cfg$seed + 622L)
  idx <- replicate(cfg$permutations, sample.int(sum(ok)))
  rp <- qr.resid(Q, (tf(x[ok]) - rx) + matrix(rx[idx], nrow = sum(ok)))
  do.call(rbind, lapply(seq_len(ncol(Y)), function(j) {
    ry <- qr.resid(Q, tf(Y[ok, j])); r <- sum(rx * ry) / sqrt(sum(rx^2) * sum(ry^2))
    nul <- abs(drop(crossprod(ry, rp))) / sqrt(sum(ry^2) * colSums(rp^2))
    data.frame(feature = colnames(Y)[j], method = method, adjusted = paste(adjust, collapse = "/"), n = sum(ok), rho = r,
               t_P = 2 * pt(-abs(r) * sqrt(df / max(1 - r^2, 1e-14)), df), P = (1 + sum(nul >= abs(r) - 1e-12)) / (cfg$permutations + 1))
  }))
}
