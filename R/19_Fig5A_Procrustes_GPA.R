## =============================================================================
##  19_Fig5A_Procrustes_GPA.R
##  ---------------------------------------------------------------------------
##  PANEL       Figure 5A  four-omics generalized Procrustes consensus (main plot)
##                         and the three pairwise Procrustes superimpositions
##                         (species-KO, KO-metabolite, metabolite-host) as insets
##
##  PURPOSE     Ask whether the sample configurations of the four data layers
##              agree: do samples that are similar in one layer tend to be
##              similar in the others?
##
##  INPUT       output/cache/four_omics_inputs.RData   (script 18)
##
##  OUTPUT      output/figures/Fig5A_procrustes_GPA.{svg,png}
##              output/tables/Fig5A_pairwise_procrustes.csv
##              output/tables/Fig5A_GPA_summary.csv, Fig5A_GPA_blocks.csv
##
##  METHODS
##   * Pairwise: symmetric Procrustes rotation (vegan::procrustes) of the first
##     <= 5 PCoA axes of the two layers on their overlapping samples;
##     correlation r = sqrt(1 - m^2).  Significance by a patient-block
##     permutation test (9,999 permutations): whole patients (with all their
##     visits) are exchanged only among patients with the same
##     visit-availability pattern, so baseline and after-RT samples of one
##     patient stay together - the same null as the GPA test below.  The
##     unrestricted PROTEST P (vegan::protest, 9,999 sample permutations) is
##     kept in the table as `P_unrestricted` for reference.  95 % interval =
##     r +/- 1.96 x SD of r over 999 patient-level bootstrap resamples (normal
##     approximation).
##   * Four-layer GPA (22 samples with all layers): configurations centred and
##     scaled to unit total variance, iteratively rotated to their mean
##     (consensus) until the residual changes < 1e-10.  Agreement = mean cosine
##     similarity between each aligned layer and the consensus; the
##     permutation P (9,999) permutes patients within identical
##     visit-availability patterns for layers 2-4 and compares the total
##     residual sum of squares.
##
##  CHANGE LOG
##   * v2 (revision): the pairwise P values use the patient-block permutation
##     instead of unrestricted PROTEST permutations, for consistency with the
##     GPA test and with the repeated-measures design.  r and the bootstrap
##     intervals are unchanged; the original unrestricted P is reported as
##     `P_unrestricted` (0.0001 / 0.0015 / 0.1403 in the submitted panel).
##
##  R PACKAGES  vegan, ggplot2, patchwork
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(vegan)
  library(patchwork)
})
report_versions(c("vegan", "ggplot2", "patchwork"))

load(cache_file("four_omics_inputs.RData"))                      # coherence (script 18)

n_protest <- 9999; n_boot <- 999; n_gpa_perm <- 9999          # <-- permutation / bootstrap sizes
layer_order  <- c("species", "ko", "metabolite", "host")
layer_labels <- c(species = "Species", ko = "KEGG ortholog", metabolite = "Metabolite", host = "Host RNA-seq")
legend_order <- c("Species", "Metabolite", "KEGG ortholog", "Host RNA-seq")            # legend order of the submitted panel
layer_colors <- c(Species = "#79B3A3", `KEGG ortholog` = "#D8A46F", Metabolite = "#82A8C7", `Host RNA-seq` = "#B29AC6")
layer_shapes <- c(Species = 21, `KEGG ortholog` = 24, Metabolite = 22, `Host RNA-seq` = 23)

## Patient-block permutation (shared by the pairwise tests and the GPA test):
## patients are exchanged only within identical visit-availability patterns
## ("Before", "Before|Ongoing", "Ongoing"), so every visit of a patient moves
## together and the number of baseline / after-RT rows stays fixed.  Returns
## the permuted row index for a sample table `meta` (SubjectID, Timepoint).
permute_rows <- function(meta) {
  pattern <- tapply(meta$Timepoint, meta$SubjectID, function(z) paste(sort(unique(z)), collapse = "|"))
  mapping <- setNames(names(pattern), names(pattern))
  for (p in unique(pattern)) { s <- names(pattern)[pattern == p]; mapping[s] <- sample(s) }
  match(paste(mapping[meta$SubjectID], meta$Timepoint), paste(meta$SubjectID, meta$Timepoint))
}
procrustes_r <- function(x, y) sqrt(max(0, 1 - procrustes(x, y, symmetric = TRUE)$ss))


## ---- 1. Pairwise Procrustes: species-KO, KO-metabolite, metabolite-host -------------------

pair_names <- c("species_ko", "ko_metabolite", "metabolite_host")
pair_seed_index <- c(1, 4, 6)                                    # seeds as in the original run
pair_stats <- list(); pair_plots <- list()
for (j in seq_along(pair_names)) {
  pair <- coherence$pairs[[pair_names[j]]]
  bx <- layer_order[j]; by <- layer_order[j + 1]
  x <- pair[[paste0(bx, "_pcoa_scores")]]; y <- pair[[paste0(by, "_pcoa_scores")]]
  ids <- pair$sample_meta$SampleID
  k <- min(5, ncol(x), ncol(y), length(ids) - 1)
  x <- x[ids, 1:k]; y <- y[ids, 1:k]
  meta_pair <- pair$sample_meta                                  # rows in the order of `ids`
  stopifnot(identical(as.character(meta_pair$SampleID), as.character(ids)))

  fit <- procrustes(x, y, symmetric = TRUE)
  r_obs <- procrustes_r(x, y)

  ## unrestricted PROTEST (reference only; the P of the submitted panel)
  set.seed(20260718 + pair_seed_index[j] * 100 + 1)
  test_free <- protest(x, y, permutations = n_protest, symmetric = TRUE)
  stopifnot(abs(test_free$t0 - r_obs) < 1e-10)

  ## patient-block permutation test (reported P): permute the rows of y,
  ## whole patients within identical visit-availability patterns
  set.seed(20260718 + pair_seed_index[j] * 100 + 2)
  perm_r <- replicate(n_protest, procrustes_r(x, y[permute_rows(meta_pair), , drop = FALSE]))
  p_block <- (1 + sum(perm_r >= r_obs)) / (1 + n_protest)

  ## patient-level bootstrap of r
  patient <- meta_pair$SubjectID
  set.seed(20267718 + pair_seed_index[j])
  boot_r <- replicate(n_boot, {
    rows <- unlist(lapply(sample(unique(patient), replace = TRUE), function(id) which(patient == id)))
    procrustes_r(x[rows, , drop = FALSE], y[rows, , drop = FALSE])
  })
  ci <- pmax(0, pmin(1, r_obs + c(-1, 1) * qnorm(0.975) * sd(boot_r)))
  pair_stats[[j]] <- data.frame(Pair = pair_names[j], N = length(ids), N_patients = length(unique(patient)), Axes = k,
                                r = r_obs, CI_low = ci[1], CI_high = ci[2],
                                P = p_block, P_unrestricted = test_free$signif, permutations = n_protest)

  pts <- rbind(data.frame(x = fit$X[, 1], y = fit$X[, 2], Layer = layer_labels[[bx]]),
               data.frame(x = fit$Yrot[, 1], y = fit$Yrot[, 2], Layer = layer_labels[[by]]))
  seg <- data.frame(x = fit$X[, 1], y = fit$X[, 2], xend = fit$Yrot[, 1], yend = fit$Yrot[, 2])
  pair_plots[[j]] <- ggplot(pts, aes(x, y)) +
    geom_segment(data = seg, aes(xend = xend, yend = yend), colour = "grey75", linewidth = 0.3) +
    geom_point(aes(fill = Layer, shape = Layer), colour = "grey30", size = 2) +
    scale_fill_manual(values = layer_colors, guide = "none") + scale_shape_manual(values = layer_shapes, guide = "none") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 2.2,
             label = sprintf("r = %.2f [95%% CI %.2f–%.2f]\np %s, n = %d", r_obs, ci[1], ci[2],
                             if (p_block < 0.001) "< 0.001" else sprintf("= %.*f", if (p_block < 0.01) 3 else 2, p_block), length(ids))) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.45))) +                  # room for the text above the points
    scale_x_continuous(expand = expansion(mult = 0.08)) +
    coord_equal() +
    labs(x = NULL, y = paste(layer_labels[[bx]], layer_labels[[by]], sep = " – ")) +
    theme_tnt(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), axis.title.y = element_text(size = 7))
}
pair_stats <- do.call(rbind, pair_stats)
print(pair_stats, row.names = FALSE, digits = 4)
save_table(pair_stats, "Fig5A_pairwise_procrustes")


## ---- 2. Generalized Procrustes analysis of the four layers ----------------------------------

## Configurations: for every layer the pair-specific PCoA scores that contain
## all common samples, preferring the largest sample set (ties -> alphabetical)
common <- coherence$four_block$sample_meta$SampleID
gpa_source <- sapply(layer_order, function(b) {
  cand <- Filter(function(p) b %in% p$block_names && all(common %in% rownames(p[[paste0(b, "_pcoa_scores")]])), coherence$pairs)
  info <- data.frame(pair = names(cand), axes = sapply(cand, function(p) ncol(p[[paste0(b, "_pcoa_scores")]])), n = sapply(cand, function(p) p$n))
  info$pair[order(-info$axes, -info$n, info$pair)][1]
})
k <- min(5, length(common) - 1, sapply(layer_order, function(b) ncol(coherence$pairs[[gpa_source[[b]]]][[paste0(b, "_pcoa_scores")]])))
configs <- lapply(layer_order, function(b) coherence$pairs[[gpa_source[[b]]]][[paste0(b, "_pcoa_scores")]][common, 1:k])
names(configs) <- layer_order

## GPA: centre + unit-norm each configuration, rotate iteratively to the consensus
fit_gpa <- function(configs, tol = 1e-10, max_iter = 500, pairwise = TRUE) {
  std <- lapply(configs, function(x) { x <- scale(x, scale = FALSE); x / sqrt(sum(x^2)) })
  normalise <- function(m) { m <- scale(m, scale = FALSE); m / sqrt(sum(m^2)) }
  rotate_to <- function(x, target) { s <- svd(crossprod(x, target)); x %*% s$u %*% t(s$v) }
  consensus <- normalise(Reduce(`+`, std) / length(std)); prev <- Inf; converged <- FALSE; iter <- max_iter
  for (i in seq_len(max_iter)) {
    aligned <- lapply(std, rotate_to, target = consensus)
    new_consensus <- normalise(Reduce(`+`, aligned) / length(aligned))
    loss <- sum(sapply(aligned, function(x) sum((x - new_consensus)^2)))
    consensus <- new_consensus
    if (is.finite(prev) && abs(prev - loss) < tol) { converged <- TRUE; iter <- i; break }
    prev <- loss
  }
  aligned <- lapply(std, rotate_to, target = consensus)
  rot <- svd(consensus, nu = 0)$v                                 # principal axes of the consensus
  consensus <- consensus %*% rot; aligned <- lapply(aligned, function(x) x %*% rot)
  residual <- sapply(aligned, function(x) sum((x - consensus)^2))
  similarity <- sapply(aligned, function(x) sum(x * consensus) / sqrt(sum(x^2) * sum(consensus^2)))
  pw <- matrix(1, length(aligned), length(aligned), dimnames = list(names(aligned), names(aligned)))
  if (pairwise) for (a in 1:(length(aligned) - 1)) for (b in (a + 1):length(aligned)) {
    pw[a, b] <- pw[b, a] <- sqrt(max(0, 1 - procrustes(aligned[[a]], aligned[[b]], symmetric = TRUE)$ss))
  }
  list(aligned = aligned, consensus = consensus, residual_ss = residual, total_residual_ss = sum(residual),
       similarity = similarity, agreement = mean(similarity), pairwise_r = pw, mean_pairwise_r = mean(pw[upper.tri(pw)]),
       axis_fraction = colSums(consensus^2) / sum(consensus^2), converged = converged, iterations = iter)
}
gpa <- fit_gpa(configs)
cat(sprintf("GPA: n = %d, axes = %d, agreement = %.4f, mean pairwise r = %.4f, converged in %d iterations\n",
            length(common), k, gpa$agreement, gpa$mean_pairwise_r, gpa$iterations))

## Permutation test: patients permuted within identical visit-availability
## patterns (permute_rows, defined above - the same scheme as the pairwise tests)
meta4 <- coherence$four_block$sample_meta
stopifnot(identical(as.character(meta4$SampleID), as.character(common)))
set.seed(20260728)
perm_ss <- replicate(n_gpa_perm, {
  pc <- configs
  for (b in 2:length(pc)) pc[[b]] <- pc[[b]][permute_rows(meta4), , drop = FALSE]
  fit_gpa(pc, pairwise = FALSE)$total_residual_ss
})
gpa_p <- (1 + sum(perm_ss <= gpa$total_residual_ss)) / (1 + length(perm_ss))
cat(sprintf("GPA permutation P = %.4f\n", gpa_p))

save_table(data.frame(n_samples = length(common), axes = k, agreement = gpa$agreement, mean_pairwise_r = gpa$mean_pairwise_r,
                      total_residual_ss = gpa$total_residual_ss, permutation_P = gpa_p, permutations = n_gpa_perm,
                      converged = gpa$converged, iterations = gpa$iterations), "Fig5A_GPA_summary")
save_table(data.frame(layer = layer_order, source_pair = unname(gpa_source[layer_order]), residual_ss = gpa$residual_ss[layer_order],
                      consensus_similarity = gpa$similarity[layer_order]), "Fig5A_GPA_blocks")


## ---- 3. Figure 5A -------------------------------------------------------------------------

consensus_df <- data.frame(SampleID = common, C1 = gpa$consensus[, 1], C2 = gpa$consensus[, 2])
points_df <- do.call(rbind, lapply(layer_order, function(b) data.frame(SampleID = common, A1 = gpa$aligned[[b]][, 1], A2 = gpa$aligned[[b]][, 2],
                                                                        Layer = layer_labels[[b]])))
points_df$Layer <- factor(points_df$Layer, levels = legend_order)
segments_df <- merge(points_df, consensus_df, by = "SampleID")
summary_text <- sprintf("Agreement = %.2f\nMean all-pair r = %.2f\nPermutation, p %s", gpa$agreement, gpa$mean_pairwise_r,
                        ifelse(gpa_p < 0.001, "< 0.001", sprintf("= %.3f", gpa_p)))

p_gpa <- ggplot() +
  geom_segment(data = segments_df, aes(C1, C2, xend = A1, yend = A2, colour = Layer), linewidth = 0.35, alpha = 0.5) +
  geom_point(data = consensus_df, aes(C1, C2), shape = 4, size = 1.4, stroke = 0.5, colour = "grey40") +
  geom_point(data = points_df, aes(A1, A2, shape = Layer, fill = Layer), size = 2.6, colour = "grey25", stroke = 0.5) +
  annotate("label", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.3, label = summary_text, size = 2.8, lineheight = 1,
           fill = "white") +
  scale_shape_manual(values = layer_shapes, name = "Data type") +
  scale_fill_manual(values = layer_colors, name = "Data type") +
  scale_colour_manual(values = layer_colors, guide = "none") +
  labs(x = sprintf("Consensus axis 1 (%.1f%%)", 100 * gpa$axis_fraction[1]), y = sprintf("Consensus axis 2 (%.1f%%)", 100 * gpa$axis_fraction[2])) +
  coord_fixed(clip = "off") +
  theme_tnt(base_size = 10, legend = "bottom") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.title = element_text(face = "bold", size = 9), legend.text = element_text(size = 8.5),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = 0.4))

p_5a <- (p_gpa | wrap_plots(pair_plots, ncol = 1)) + plot_layout(widths = c(3, 1.15))
save_panel(p_5a, "Fig5A_procrustes_GPA", width = 8.2, height = 6.2)

message("Done: 19_Fig5A_Procrustes_GPA.R")
