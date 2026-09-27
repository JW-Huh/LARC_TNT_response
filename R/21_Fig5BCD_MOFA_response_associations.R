## =============================================================================
##  21_Fig5BCD_MOFA_response_associations.R
##  ---------------------------------------------------------------------------
##  PANELS      Figure 5B  every MOFA factor: Hedges' g (pCR - non-pCR) overall /
##                         baseline / after RT, variance explained per view, total
##              Figure 5C  Factor 7 scores by response
##              Figure 5D  Factor 7 x Factor 11 plane with PERMANOVA
##
##  PURPOSE     Test which latent factors of the joint four-omics model separate
##              pCR from non-pCR patients, and characterise the leading factor.
##
##  INPUT       output/cache/mofa_model.RData   (script 20)
##
##  OUTPUT      output/figures/Fig5B_MOFA_factor_associations, Fig5C_Factor7_scores,
##              Fig5D_Factor7_Factor11  (.svg/.png)
##              output/tables/Fig5B_factor_associations.csv, Fig5B_period_effects.csv,
##              Fig5B_variance_explained.csv, Fig5D_factor_pair_tests.csv
##              output/cache/mofa_associations.RData
##
##  METHODS     Factor scores are standardised (all 42 observations).  For each
##              factor: random-intercept model score ~ Response (nlme::lme,
##              patient random effect); Hedges' g on patient means (bootstrap
##              95 % CI, 999 resamples); permutation P from the GLS F statistic
##              with 999 permutations of the patient-level response label
##              restricted to patients with the same visit pattern (baseline
##              only / both visits); BH q over the 14 factors and a max-F
##              family-wise P.  Variance explained is reconstructed from scores
##              x weights on the model-scaled data and checked against the
##              values stored by MOFA2.  Factor pair: PERMANOVA-type F on
##              Euclidean distance in the two-factor plane, same permutations.
##              NOTE  Factor 7: nominal permutation P = 0.011, BH q = 0.15,
##              max-F P = 0.16 - an exploratory association, not FDR-significant.
##
##  R PACKAGES  nlme, ggplot2, ggbeeswarm, patchwork
## =============================================================================

source("R/_common.R")
source("R/_helpers_mofa.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggbeeswarm)
  library(patchwork)
})
report_versions(c("nlme", "dplyr", "ggplot2", "ggbeeswarm", "patchwork"))

load(cache_file("mofa_model.RData"))                                  # mofa_input, mofa_model


## ---- 1. Associations --------------------------------------------------------------------

a <- analyze_association(mofa_model$extracted, mofa_input$metadata)
effects <- period_effects(a)
print(a$association[order(a$association$permutation_P), c("factor", "beta", "g", "permutation_P", "permutation_q", "maxF_P")], row.names = FALSE, digits = 3)
cat("Display factors:", paste(a$axes, collapse = ", "), "\n")
save_table(a$association, "Fig5B_factor_associations")
save_table(effects, "Fig5B_period_effects")
save_table(a$pair_results, "Fig5D_factor_pair_tests")

## variance explained: reconstructed and checked against the MOFA2 values
ve <- variance_explained(mofa_model$extracted, mofa_input$values)
stored <- mofa_model$extracted$variance
chk <- merge(stored[, c("view", "factor", "r2")], ve$components[, c("view", "factor", "r2")], by = c("view", "factor"), suffixes = c("_mofa", "_rebuilt"))
stopifnot(max(abs(chk$r2_mofa - chk$r2_rebuilt)) < 5e-4)
total_r2 <- ve$total$r2[ve$total$factor == "ALL_FACTORS"]
cat(sprintf("Joint model R^2 = %.3f (all factors, all views)\n", total_r2))
save_table(ve$components, "Fig5B_variance_explained")
save(a, effects, ve, file = cache_file("mofa_associations.RData"))


## ---- 2. Figure 5B ------------------------------------------------------------------------

factor_order <- a$association$factor[order(a$association$beta)]     # strongest pCR association on top
factor_label <- function(x) sub("Factor", "Factor ", x)
effects <- effects %>% dplyr::mutate(factor = factor(factor, levels = factor_order),
                                     period = factor(period, levels = c("Overall", "Baseline", "After RT")),
                                     tier = factor(forest_tier(P), levels = c("P < 0.05", "0.05 <= P < 0.10", "NS")),
                                     y = as.numeric(factor) + c(Overall = 0.25, Baseline = 0, `After RT` = -0.25)[as.character(period)])
tier_colors <- c("P < 0.05" = "#258F78", "0.05 <= P < 0.10" = "#8ED0BC", "NS" = "#6E6E6E")

p_forest <- ggplot(effects, aes(g, y, colour = tier, shape = period)) +
  geom_hline(yintercept = seq_along(factor_order) - 0.5, colour = "grey92", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2, linewidth = 0.3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, linewidth = 0.45, orientation = "y") +
  geom_point(aes(size = period == "Overall"), fill = "white", stroke = 0.6) +
  scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 2.2), guide = "none") +
  scale_colour_manual(values = tier_colors, name = "Sig.", labels = c(expression(italic(P) < 0.05), expression(0.05 <= italic(P) * " < 0.1"), "NS")) +
  scale_shape_manual(values = c(Overall = 21, Baseline = 22, `After RT` = 23), name = "Time") +
  scale_y_continuous(breaks = seq_along(factor_order), labels = factor_label(factor_order), expand = expansion(add = 0.6)) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(x = expression("Hedges' " * italic(g)), y = NULL) +
  theme_tnt(base_size = 10, legend = "bottom") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.box = "vertical",
        legend.background = element_rect(colour = "black", linewidth = 0.3), legend.margin = margin(2, 4, 2, 4))

view_labels <- c(species = "Species", metabolite = "Metabolite", ko = "KO", host = "Host RNA")
v <- stored %>% dplyr::filter(factor %in% factor_order) %>%
  dplyr::mutate(factor = factor(factor, levels = factor_order), view = factor(view_labels[view], levels = view_labels), pct = 100 * r2)
p_tiles <- ggplot(v, aes(view, factor, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f", pct), colour = pct > 0.52 * max(pct)), size = 3) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradient(low = "white", high = "#142D4E", name = "Within-view variance explained (%)",
                      guide = guide_colourbar(barwidth = unit(4, "cm"), barheight = unit(0.3, "cm"), title.position = "top")) +
  scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = expansion(add = 0.6)) +
  labs(x = NULL, y = NULL) +
  theme_tnt(base_size = 10, legend = "bottom") +
  theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_rect(fill = "grey95", colour = NA))

t <- ve$total %>% dplyr::filter(factor != "ALL_FACTORS") %>% dplyr::mutate(factor = factor(factor, levels = factor_order))
p_total <- ggplot(t, aes(100 * r2, factor)) +
  geom_col(fill = "grey20", width = 0.7) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3.2, label = sprintf("Joint model\nR² = %.3f", total_r2)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) + scale_y_discrete(expand = expansion(add = 0.6)) +
  labs(x = "Total variance explained (%)", y = NULL) +
  theme_tnt(base_size = 10) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

save_panel(p_forest + p_tiles + p_total + plot_layout(widths = c(2.2, 1.5, 1.2)), "Fig5B_MOFA_factor_associations", width = 9.5, height = 6.2)


## ---- 3. Figure 5C / 5D: the leading factors ------------------------------------------------

d <- a$metadata
d$Response_label <- factor(ifelse(d$Response == 1, "pCR", "non-pCR"), levels = response_levels)
d$Visit <- factor(d$Timepoint, levels = c("Before", "Ongoing"), labels = c("Baseline", "After RT"))
d$X <- a$Z[, a$axes[1]] * factor_sign(a, a$axes[1])                 # oriented so that pCR scores higher
d$Y <- a$Z[, a$axes[2]] * factor_sign(a, a$axes[2])
p_lead <- a$association$permutation_P[a$association$factor == a$axes[1]]

y_top <- max(d$X) + 0.15 * diff(range(d$X))
p_5c <- ggplot(d, aes(Response_label, X, colour = Response_label)) +
  geom_boxplot(width = 0.45, fill = NA, outlier.shape = NA, coef = 0, linewidth = 0.6) +
  geom_quasirandom(width = 0.15, size = 2.2) +
  annotate("segment", x = 1, xend = 2, y = y_top, yend = y_top, linewidth = 0.4) +
  annotate("text", x = 1.5, y = y_top, vjust = -0.6, size = 3.5, label = sprintf("italic(P) == %.3f", p_lead), parse = TRUE) +
  scale_colour_manual(values = response_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "Factor score", title = factor_label(a$axes[1])) +
  theme_tnt(base_size = 11, legend = "none") +
  theme(axis.text.x = element_text(size = rel(1.1)))
save_panel(p_5c, "Fig5C_Factor7_scores", width = 2.6, height = 4)

pair <- a$pair_results[a$best_pair, ]
hulls <- do.call(rbind, lapply(split(d, d$Response_label), function(z) z[chull(z$X, z$Y), ]))
centroids <- d %>% dplyr::group_by(Response_label) %>% dplyr::summarise(X = mean(X), Y = mean(Y), .groups = "drop")
p_5d <- ggplot(d, aes(X, Y)) +
  geom_polygon(data = hulls, aes(fill = Response_label, colour = Response_label), alpha = 0.15, linewidth = 0.7, show.legend = FALSE) +
  geom_point(aes(colour = Response_label, shape = Visit), size = 2) +
  geom_point(data = centroids, aes(fill = Response_label, shape = "Centroid"), size = 4.5, colour = "grey20") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, size = 3.2, parse = TRUE,
           label = sprintf('atop("PERMANOVA", paste(R^2 == %.3f, "; ", italic(P) == %.3f))', pair$R2, pair$P)) +
  scale_colour_manual(values = response_colors, name = NULL) +
  scale_fill_manual(values = response_colors, guide = "none") +
  scale_shape_manual(values = c(Baseline = 16, `After RT` = 15, Centroid = 23), breaks = c("Baseline", "After RT", "Centroid"), name = NULL) +
  guides(shape = guide_legend(order = 1, override.aes = list(fill = "grey70", colour = "grey30", size = c(2, 2, 3))),
         colour = guide_legend(order = 2, override.aes = list(shape = 16, size = 2))) +
  labs(x = paste(factor_label(a$axes[1]), "score"), y = paste(factor_label(a$axes[2]), "score")) +
  theme_tnt(base_size = 10, legend = "inside") +
  theme(legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0), legend.spacing.y = unit(0, "cm"),
        legend.key.spacing.y = unit(0, "pt"), legend.background = element_blank())
save_panel(p_5d, "Fig5D_Factor7_Factor11", width = 4.2, height = 4)

message("Done: 21_Fig5BCD_MOFA_response_associations.R")
