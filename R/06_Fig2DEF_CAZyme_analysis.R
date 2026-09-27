## =============================================================================
##  06_Fig2DEF_CAZyme_analysis.R
##  ---------------------------------------------------------------------------
##  PANELS      Figure 2D  CAZyme class-level shift (median pairwise log2 difference
##                         of the subfamilies in each class, with bootstrap interval)
##              Figure 2E  volcano plot of CAZyme subfamilies
##              Figure 2F  starch/alpha-glucan capacity score and five
##                         response-associated subfamilies (box + points)
##
##  PURPOSE     Test whether carbohydrate-active enzyme (CAZyme) gene content of
##              the baseline microbiome differs between pCR and non-pCR.
##
##  INPUT       input/microbiome/cazy_subfamily_tpm.tsv   CAZy subfamily TPM
##              input/microbiome/cazy_family_tpm.tsv      CAZy family TPM
##                  (dbCAN annotation of assembled genes; transcripts per million)
##              input/metadata/stool_samples.csv, patients.csv
##
##  OUTPUT      output/figures/Fig2D_CAZyme_class_shift, Fig2E_CAZyme_volcano,
##              Fig2F_CAZyme_boxplots  (.svg/.png)
##              output/tables/Fig2E_CAZyme_subfamily_statistics.csv,
##              Fig2D_CAZyme_class_summary.csv, Fig2F_CAZyme_panel_tests.csv
##
##  METHODS
##   * Subfamily filter: prevalence >= 30 % of baseline samples and mean TPM >= 1.
##   * Effect size: median of all pairwise log2(TPM + 1) differences between the
##     11 pCR and 15 non-pCR samples (a Hodges-Lehmann-type location shift,
##     positive = higher in pCR).  P: two-sided Wilcoxon rank-sum test (normal
##     approximation, exact = FALSE).
##   * Class shift (2D): median of the subfamily effects within each CAZy class
##     (GH, GT, PL, CE, CBM, AA); the interval is the 2.5-97.5 % range of the
##     median over 2,000 bootstrap resamples of the SUBFAMILIES (not of
##     patients) - a descriptive interval of the feature-effect distribution.
##   * Starch score (2F): family-level log10(TPM + 1) of nine starch/alpha-glucan
##     families (GH13, GH65, GH77, GH97, CBM20, CBM25, CBM26, CBM34, CBM48),
##     z-scored per family and averaged per sample.  Box-plot panels use the
##     default (exact when possible) two-sided Wilcoxon test.
##
##  R PACKAGES  ggplot2, ggrepel, ggbeeswarm, patchwork, dplyr
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggrepel)
  library(ggbeeswarm)
  library(patchwork)
})
report_versions(c("ggplot2", "ggrepel", "ggbeeswarm", "patchwork"))


## ---- 1. Baseline CAZyme subfamily table ----------------------------------------

baseline <- read_stool_samples(timepoint = "Before")
subfam <- read_feature_table("cazy_subfamily_tpm.tsv", "CAZy_subfamily", samples = baseline$SampleID)
stopifnot(all(is.finite(subfam)), all(subfam >= 0))
is_pcr <- baseline$Response == "pCR"

prevalence_cut <- 0.30                                    # <-- subfamily filters
mean_tpm_cut   <- 1

log_sub <- log2(subfam + 1)
stats <- data.frame(
  Feature    = rownames(subfam),
  Class      = sub("[0-9].*", "", rownames(subfam)),     # "GH13_e104" -> "GH"
  Prevalence = rowMeans(subfam > 0),
  Mean_TPM   = rowMeans(subfam),
  Effect     = apply(log_sub, 1, function(v) median(as.vector(outer(v[is_pcr], v[!is_pcr], "-")))),
  row.names  = NULL
)
stats <- stats[stats$Prevalence >= prevalence_cut & stats$Mean_TPM >= mean_tpm_cut, ]
stats$P <- apply(log_sub[stats$Feature, ], 1, function(v) wilcox_p(v, baseline$Response, exact = FALSE))
stats$Class <- factor(stats$Class, levels = c("GH", "GT", "PL", "CE", "CBM", "AA"))
stats <- stats[order(stats$P, stats$Feature), ]
cat("Subfamilies retained:", nrow(stats), "; P < 0.05:", sum(stats$P < 0.05), "\n")
save_table(stats, "Fig2E_CAZyme_subfamily_statistics")


## ---- 2. Figure 2D: class-level shift ---------------------------------------------

set.seed(123)                                             # <-- bootstrap seed
n_boot <- 2000
class_summary <- do.call(rbind, lapply(levels(stats$Class), function(cl) {
  eff <- stats$Effect[stats$Class == cl]
  boot <- replicate(n_boot, median(sample(eff, replace = TRUE)))
  data.frame(Class = cl, N = length(eff), Effect = median(eff),
             Lower = unname(quantile(boot, 0.025)), Upper = unname(quantile(boot, 0.975)))
}))
class_summary <- class_summary[order(class_summary$Effect), ]
class_summary$Class <- factor(class_summary$Class, levels = class_summary$Class)   # largest shift on top
print(class_summary, row.names = FALSE, digits = 3)
save_table(class_summary, "Fig2D_CAZyme_class_summary")

## Design (as submitted): open circles with black interval lines, pink (non-pCR)
## and teal (pCR) half-panels, "n = " labels right of each interval.
x_lim <- c(-0.23, 0.23)
p_2d <- ggplot(class_summary, aes(Effect, Class)) +
  annotate("rect", xmin = x_lim[1], xmax = 0, ymin = -Inf, ymax = Inf, fill = response_colors[["non-pCR"]], alpha = 0.12) +
  annotate("rect", xmin = 0, xmax = x_lim[2], ymin = -Inf, ymax = Inf, fill = response_colors[["pCR"]], alpha = 0.12) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_segment(aes(x = Lower, xend = Upper, yend = Class), linewidth = 0.9) +
  geom_point(shape = 21, fill = "white", size = 3.6, stroke = 1) +
  geom_text(aes(x = Upper, label = paste0("n = ", N)), hjust = -0.25, size = 3.2) +
  annotate("text", x = x_lim[1] + 0.005, y = Inf, label = "Higher in non-pCR", hjust = 0, vjust = 1.8,
           size = 2.9, colour = response_colors[["non-pCR"]]) +
  annotate("text", x = x_lim[2] - 0.005, y = Inf, label = "Higher in pCR", hjust = 1, vjust = 1.8,
           size = 2.9, colour = response_colors[["pCR"]]) +
  scale_x_continuous(limits = x_lim, breaks = seq(-0.2, 0.2, 0.1), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 1.2))) +
  labs(x = expression(atop("CAZyme subfamily shift", "(Pairwise " * log[2] * " fold-change)")), y = NULL) +
  theme_tnt(base_size = 11, legend = "none") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_panel(p_2d, "Fig2D_CAZyme_class_shift", width = 4.2, height = 4.2)


## ---- 3. Figure 2E: subfamily volcano ----------------------------------------------

stats$Direction <- ifelse(stats$P < 0.05, ifelse(stats$Effect > 0, "pCR-enriched", "non-pCR-enriched"), "ns")

## Subfamilies labelled in the submitted panel  <-- edit to label other features
label_features <- c("AA4", "GH95_e26", "GH78_e165", "GH2_e57", "CE2_e37", "CBM50_e1290",
                    "GH13_e104", "CE4_e125", "GH1_e175", "CBM34_e9", "GH13_e374")

p_2e <- ggplot(stats, aes(Effect, -log10(P))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(data = filter(stats, Direction == "ns"), colour = "grey78", size = 1.1, shape = 18) +
  geom_point(data = filter(stats, Direction != "ns"), aes(colour = Direction), size = 2.2) +
  geom_text_repel(data = filter(stats, Feature %in% label_features), aes(label = Feature, colour = Direction),
                  size = 2.8, seed = 123, max.overlaps = Inf, show.legend = FALSE, min.segment.length = 0.2) +
  scale_colour_manual(values = c("pCR-enriched" = response_colors[["pCR"]],
                                 "non-pCR-enriched" = response_colors[["non-pCR"]]),
                      breaks = c("non-pCR-enriched", "pCR-enriched"), name = NULL) +
  labs(x = expression("Pairwise " * log[2] * " fold-change"), y = expression(-Log[10] * " (" * italic(p) * ")")) +
  theme_tnt(base_size = 11, legend = "inside") +
  theme(legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_blank())
save_panel(p_2e, "Fig2E_CAZyme_volcano", width = 4.2, height = 4.2)


## ---- 4. Figure 2F: starch/alpha-glucan score and five subfamilies ------------------

family <- read_feature_table("cazy_family_tpm.tsv", "CAZy_family", samples = baseline$SampleID)
starch_families <- c("GH13", "GH65", "GH77", "GH97", "CBM20", "CBM25", "CBM26", "CBM34", "CBM48")   # <-- score members
stopifnot(all(starch_families %in% rownames(family)))

starch_z <- t(scale(t(log10(family[starch_families, ] + 1))))   # z-score each family across samples
starch_z[!is.finite(starch_z)] <- 0                             # a constant family contributes 0
starch_score <- colMeans(starch_z)

## Representative subfamilies with a one-line functional description
representatives <- c(
  "GH2_e57"     = "(GH2-family glycoside hydrolase)",
  "CE2_e37"     = "(CE2-family carbohydrate esterase)",
  "CBM50_e1290" = "(Peptidoglycan- or chitin-binding)",
  "AA4"         = "(Vanillyl-alcohol oxidase)",
  "GH95_e26"    = "(α-L-fucosidase-related activity)"
)
stopifnot(all(names(representatives) %in% rownames(subfam)))

panel_data <- list("Starch/α-glucan score" = list(value = starch_score, unit = "Capacity score", note = ""))
for (f in names(representatives)) {
  panel_data[[f]] <- list(value = log10(subfam[f, ] + 1), unit = expression(log[10] * "(TPM + 1)"), note = representatives[[f]])
}

box_panel <- function(name, d) {
  df <- data.frame(Response = baseline$Response, Value = d$value)
  p <- wilcox_p(df$Value, df$Response)                            # default exact-when-possible Wilcoxon
  ggplot(df, aes(Response, Value, colour = Response)) +
    geom_boxplot(width = 0.35, fill = NA, outlier.shape = NA, coef = 0, linewidth = 0.6) +
    geom_quasirandom(width = 0.15, size = 1.8) +
    annotate("text", x = 1.5, y = Inf, vjust = 1.5, size = 3.1,
             label = paste0("italic(p) == ", formatC(p, format = "g", digits = 2)), parse = TRUE) +
    scale_colour_manual(values = response_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = NULL, y = d$unit, title = name, subtitle = d$note) +
    theme_tnt(base_size = 9, legend = "none") +
    theme(plot.subtitle = element_text(hjust = 0.5, size = 6, colour = "grey30"),
          plot.title = element_text(size = 9.5))
}
panels <- Map(box_panel, names(panel_data), panel_data)
panel_tests <- data.frame(Panel = names(panel_data),
                          P = sapply(panel_data, function(d) wilcox_p(d$value, baseline$Response)))
print(panel_tests, row.names = FALSE)
save_table(panel_tests, "Fig2F_CAZyme_panel_tests")

save_panel(wrap_plots(panels, ncol = 3), "Fig2F_CAZyme_boxplots", width = 6.4, height = 5.2)

message("Done: 06_Fig2DEF_CAZyme_analysis.R")
