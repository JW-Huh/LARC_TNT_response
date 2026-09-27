## =============================================================================
##  07_Fig3BF_SFig3A_Metabolite_effects.R
##  ---------------------------------------------------------------------------
##  PANELS      Figure 3B   forest plot: standardised effect (Hedges' g) of the
##                          seven baseline fecal metabolites with Wilcoxon P < 0.1
##              Suppl. Fig. 3A   the same forest plot for all 35 testable metabolites
##              Figure 3F   ILA / IPA / IAA ternary display + concentration-ratio box plot
##
##  PURPOSE     Quantify baseline differences in fecal metabolite concentrations
##              (targeted LC-MS/GC-MS panel of 40 analytes) between pCR and
##              non-pCR, with an effect size and a resampling interval rather
##              than only a P value, and display the balance of the reductive
##              (ILA, IPA) versus decarboxylative (IAA) tryptophan branch.
##
##  INPUT       input/metabolome/fecal_metabolites.csv   40 metabolites (umol/g),
##                    32 samples (20 baseline / 12 after RT)
##              input/metadata/patients.csv
##
##  OUTPUT      output/figures/Fig3B_metabolite_forest, SFig3A_metabolite_forest_all,
##              Fig3F_tryptophan_ternary  (.svg/.png)
##              output/tables/Fig3B_SFig3A_metabolite_effects.csv
##              output/tables/Fig3F_tryptophan_ratio.csv
##
##  METHODS
##   * Baseline metabolome: 20 patients (8 pCR, 12 non-pCR).  Metabolites
##     detected (> 0) in >= 20 % of samples (>= 4) and with a defined Wilcoxon
##     P (non-constant) are testable: 35 of 40.
##   * P: two-sided Wilcoxon rank-sum test, normal approximation with continuity
##     correction (exact = FALSE).  Metabolites with P < 0.1 are shown in 3B.
##   * Effect size: Hedges' g = (mean_pCR - mean_nonpCR) / pooled SD, multiplied
##     by the small-sample correction 1 - 3/(4 df - 1); positive = higher in pCR.
##     95 % interval = 2.5-97.5 percentiles of g over 2,000 group-stratified
##     bootstrap resamples (seed 123).  This is a percentile bootstrap of the
##     effect size, not an inversion of the Wilcoxon test.
##   * Ternary (3F): ILA, IPA and IAA concentrations are z-scored per
##     metabolite and passed through a row-wise softmax, so each patient gets
##     three weights summing to 1 - a DISPLAY transformation; the coordinates
##     are not concentration fractions.  Hulls = convex hull per group,
##     diamonds = group centroids.  The ratio test uses raw concentrations:
##     log10((IPA + ILA + 1) / (IAA + 1)), Wilcoxon rank-sum (exact = FALSE).
##
##  R PACKAGES  dplyr, tidyr, ggplot2, patchwork
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(patchwork)
})
report_versions(c("ggplot2", "dplyr", "patchwork"))


## ---- 1. Baseline metabolite table --------------------------------------------------

metab <- read_metabolites(timepoint = "Before")                 # 20 patients x 40 metabolites
metab <- metab[order(metab$EnrolmentOrder), ]                   # patient order of the original run (affects the bootstrap draws only)
metabolites <- metabolite_columns()
is_pcr <- metab$Response == "pCR"
cat("Baseline metabolomes:", nrow(metab), "(pCR", sum(is_pcr), "/ non-pCR", sum(!is_pcr), ")\n")

prevalence_cut <- 0.20                                         # <-- detection filter
detected <- colSums(metab[, metabolites] > 0) >= ceiling(prevalence_cut * nrow(metab))
testable <- metabolites[detected]


## ---- 2. Wilcoxon P, Hedges' g and bootstrap interval ---------------------------------

n_boot <- 2000                                                  # <-- bootstrap resamples

effect_of <- function(x, y) {
  data.frame(
    n_pCR = length(x), n_nonpCR = length(y),
    Median_pCR = median(x), Median_nonpCR = median(y),
    Mean_pCR = mean(x), Mean_nonpCR = mean(y),
    P = wilcox_p(c(x, y), rep(response_levels, c(length(x), length(y))), exact = FALSE),
    Hedges_g = hedges_g(x, y)
  )
}

effects <- do.call(rbind, lapply(testable, function(m) {
  cbind(Metabolite = m, effect_of(metab[is_pcr, m], metab[!is_pcr, m]))
}))
effects <- effects[is.finite(effects$P), ]                      # constant metabolites have no P
effects <- effects[order(effects$P), ]
cat("Testable metabolites:", nrow(effects), "\n")

## Bootstrap in the same (P-ordered) sequence and with the same seed as the
## submitted analysis so that the intervals are reproduced exactly.
set.seed(123)                                                   # <-- bootstrap seed
ci <- t(sapply(effects$Metabolite, function(m) {
  x <- metab[is_pcr, m]; y <- metab[!is_pcr, m]
  g <- replicate(n_boot, hedges_g(sample(x, replace = TRUE), sample(y, replace = TRUE)))
  quantile(g, c(0.025, 0.975), na.rm = TRUE)
}))
effects$CI_low  <- ci[, 1]
effects$CI_high <- ci[, 2]
effects$Direction <- ifelse(effects$Hedges_g >= 0, "Higher in pCR", "Higher in non-pCR")
effects$Label <- gsub("_", " ", effects$Metabolite)
effects$Label[effects$Metabolite == "Indolepropionic_acid"] <- "Indole propionic acid"
effects$Label[effects$Metabolite == "Tauroursodeoxycholic_acid_Taurohyodeoxycholic_acid"] <- "TUDCA / THDCA"
effects$Label[effects$Metabolite == "gamma_Aminobutyric_acid"] <- "γ-Aminobutyric acid"
save_table(effects, "Fig3B_SFig3A_metabolite_effects")

screen_p_cut <- 0.10                                            # <-- Figure 3B inclusion threshold
selected <- effects[effects$P < screen_p_cut, ]
print(selected[, c("Metabolite", "P", "Hedges_g", "CI_low", "CI_high")], row.names = FALSE, digits = 3)


## ---- 3. Forest plot (Figure 3B, Supplementary Figure 3A) ----------------------------
## Design (as submitted): metabolite names on the y axis; group medians as
## coloured text columns left of the axis; point + interval coloured by
## direction; "g [95 % CI]" text on the right; arrows in the header.

forest_plot <- function(d, base_size = 10) {
  d <- d[order(d$Hedges_g), ]                                   # largest g at the top
  d$Label <- factor(d$Label, levels = d$Label)
  n <- nrow(d)
  x_text <- c(pCR = -4.6, nonpCR = -3.6, ci = 4.3)              # x positions of the text columns
  header_y <- n + 0.9
  ggplot(d, aes(Hedges_g, Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high, colour = Direction), width = 0.18, linewidth = 0.55, orientation = "y") +
    geom_point(aes(colour = Direction), size = 2.8) +
    geom_text(aes(x = x_text[["pCR"]], label = sprintf("%.1e", Median_pCR)), colour = response_colors[["pCR"]], size = base_size * 0.26) +
    geom_text(aes(x = x_text[["nonpCR"]], label = sprintf("%.1e", Median_nonpCR)), colour = response_colors[["non-pCR"]], size = base_size * 0.26) +
    geom_text(aes(x = x_text[["ci"]], label = sprintf("%.2f [%.2f, %.2f]", Hedges_g, CI_low, CI_high)), size = base_size * 0.26) +
    ## header row
    annotate("text", x = x_text[["pCR"]], y = header_y, label = "pCR", colour = response_colors[["pCR"]], size = base_size * 0.28) +
    annotate("text", x = x_text[["nonpCR"]], y = header_y, label = "non-pCR", colour = response_colors[["non-pCR"]], size = base_size * 0.28) +
    annotate("text", x = x_text[["ci"]], y = header_y, label = "Hedges' g [95% CI]", size = base_size * 0.28) +
    annotate("segment", x = -0.4, xend = -2.3, y = header_y, yend = header_y, colour = response_colors[["non-pCR"]],
             arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.5) +
    annotate("segment", x = 0.9, xend = 2.3, y = header_y, yend = header_y, colour = response_colors[["pCR"]],
             arrow = arrow(length = unit(0.18, "cm")), linewidth = 0.5) +
    annotate("text", x = -0.15, y = header_y, label = "non-pCR", hjust = 1, fontface = "bold", colour = response_colors[["non-pCR"]], size = base_size * 0.3) +
    annotate("text", x = 0.15, y = header_y, label = "pCR", hjust = 0, fontface = "bold", colour = response_colors[["pCR"]], size = base_size * 0.3) +
    scale_colour_manual(values = c("Higher in pCR" = response_colors[["pCR"]], "Higher in non-pCR" = response_colors[["non-pCR"]]), guide = "none") +
    scale_x_continuous(limits = c(-5.2, 5.6), breaks = c(-2.5, 0, 2.5), expand = c(0, 0)) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1.5))) +
    labs(x = "Standardized effect size", y = NULL) +
    theme_tnt(base_size = base_size, legend = "none") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size = rel(1.1)))
}
save_panel(forest_plot(selected), "Fig3B_metabolite_forest", width = 6.2, height = 4.2)
save_panel(forest_plot(effects, base_size = 9), "SFig3A_metabolite_forest_all", width = 6.6, height = 11)


## ---- 4. Figure 3F: ILA / IPA / IAA ternary display and ratio test -------------------

trp <- data.frame(Response = metab$Response,
                  ILA = metab$Indole_lactic_acid, IPA = metab$Indolepropionic_acid, IAA = metab$Indole_acetic_acid)

## Display transformation: z-score per metabolite, then row-wise softmax -> weights summing to 1
z_trp <- scale(as.matrix(trp[, c("ILA", "IPA", "IAA")]))
w <- exp(z_trp - apply(z_trp, 1, max)); w <- w / rowSums(w)
## ternary -> Cartesian (IPA bottom-left, IAA bottom-right, ILA top)
trp$x <- 0.5 * w[, "ILA"] + w[, "IAA"]
trp$y <- sqrt(3) / 2 * w[, "ILA"]

centroids <- trp %>% group_by(Response) %>% summarise(x = mean(x), y = mean(y), .groups = "drop")
hulls <- do.call(rbind, lapply(split(trp, trp$Response), function(d) d[chull(d$x, d$y), ]))
triangle <- data.frame(x = c(0, 0.5, 1, 0), y = c(0, sqrt(3) / 2, 0, 0))

## dashed grid at 20 % steps (three families of lines)
grid_lines <- do.call(rbind, lapply(seq(0.2, 0.8, 0.2), function(f) rbind(
  data.frame(x = f, y = 0, xend = 0.5 + f / 2, yend = sqrt(3) / 2 * (1 - f)),            # constant IAA... etc.
  data.frame(x = 1 - f, y = 0, xend = (1 - f) / 2, yend = sqrt(3) / 2 * (1 - f)),
  data.frame(x = f / 2, y = sqrt(3) / 2 * f, xend = 1 - f / 2, yend = sqrt(3) / 2 * f)
)))

## The submitted panel uses a slightly deeper pink for non-pCR than the rest of the paper
ternary_colors <- c("pCR" = "#5BB295", "non-pCR" = "#C65372")

p_ternary <- ggplot(trp, aes(x, y)) +
  geom_segment(data = grid_lines, aes(x, y, xend = xend, yend = yend), colour = "grey75", linetype = "dashed", linewidth = 0.3) +
  geom_path(data = triangle, linewidth = 0.7) +
  geom_polygon(data = hulls, aes(fill = Response, colour = Response), alpha = 0.15, linewidth = 0.6) +
  geom_point(aes(colour = Response), size = 2.6) +
  geom_point(data = centroids, aes(fill = Response), shape = 23, size = 5, colour = "black") +
  annotate("text", x = c(0, 0.5, 1), y = c(-0.06, 0.92, -0.06), label = c("IPA", "ILA", "IAA"), size = 5) +
  scale_fill_manual(values = ternary_colors, name = NULL) +
  scale_colour_manual(values = ternary_colors, name = NULL) +
  coord_equal(clip = "off") +
  theme_void() +
  theme(legend.position = "inside", legend.position.inside = c(0.08, 0.72), legend.text = element_text(size = 12),
        plot.margin = margin(10, 10, 10, 10))

trp$Ratio <- log10((trp$IPA + trp$ILA + 1) / (trp$IAA + 1))      # pseudocount 1 in umol/g
ratio_p <- wilcox_p(trp$Ratio, trp$Response, exact = FALSE)
cat("ILA + IPA vs IAA ratio: Wilcoxon P =", signif(ratio_p, 3), "\n")
save_table(data.frame(metab[, c("SubjectID", "SampleID", "Response")], trp[, c("ILA", "IPA", "IAA", "Ratio")]), "Fig3F_tryptophan_ratio")

p_ratio <- ggplot(trp, aes(Response, Ratio, fill = Response, colour = Response)) +
  geom_boxplot(width = 0.6, alpha = 0.25, outlier.shape = NA, linewidth = 0.6) +
  annotate("segment", x = 1, xend = 2, y = max(trp$Ratio) + 0.1, yend = max(trp$Ratio) + 0.1) +
  annotate("text", x = 1.5, y = max(trp$Ratio) + 0.1, vjust = -0.5, size = 3.4,
           label = paste0("italic(P) == ", formatC(ratio_p, format = "f", digits = 3)), parse = TRUE) +
  scale_fill_manual(values = ternary_colors, guide = "none") +
  scale_colour_manual(values = ternary_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.25))) +
  labs(x = NULL, y = NULL, title = expression(Log[10] * " [(IPA+ILA+1)/(IAA+1)]")) +
  theme_tnt(base_size = 10, legend = "none") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        plot.title = element_text(size = 9), axis.text.x = element_text(size = rel(1.1)))

p_3f <- p_ternary + inset_element(p_ratio, left = 0.68, bottom = 0.45, right = 1.0, top = 0.95)
save_panel(p_3f, "Fig3F_tryptophan_ternary", width = 6.2, height = 4.6)

message("Done: 07_Fig3BF_SFig3A_Metabolite_effects.R")
