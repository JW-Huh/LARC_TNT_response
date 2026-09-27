## =============================================================================
##  04_Fig1H_SGB_PLS_DA.R
##  ---------------------------------------------------------------------------
##  PANEL       Figure 1H  PLS-DA of baseline SGB profiles: VIP scores of the ten
##                         most discriminating SGBs per response group, with the
##                         two-component score plot as an inset
##
##  PURPOSE     Show which species-level genome bins (SGBs) drive the
##              multivariate separation of pCR and non-pCR baseline microbiomes.
##              This is a supervised, descriptive projection - it is NOT a
##              cross-validated classifier (the random-forest models of
##              scripts 09a-09c are).
##
##  INPUT       input/microbiome/metaphlan_sgb.csv    MetaPhlAn 4 SGB relative
##                                                    abundance (%), 42 samples
##              input/metadata/stool_samples.csv, input/metadata/patients.csv
##
##  OUTPUT      output/figures/Fig1H_SGB_PLS_DA_VIP.{svg,png}
##              output/tables/Fig1H_PLS_DA_VIP_scores.csv   (all SGBs)
##              output/tables/Fig1H_PLS_DA_scores.csv       (sample coordinates)
##
##  METHODS     * PLS-DA (mixOmics 6.26.0, plsda()): partial least squares
##                regression of the dummy-coded response on the SGB matrix.
##                Two components; untransformed relative abundances (%), every
##                SGB centred and scaled to unit variance (scale = TRUE, the
##                mixOmics default) so that abundant and rare SGBs contribute
##                on the same footing.
##              * VIP (variable importance in projection, mixOmics::vip()):
##                for each SGB the weighted sum of squared loadings over the
##                components, weighted by the Y-variance each component
##                explains; VIP > 1 is the usual "important" threshold.  Here
##                the maximum VIP over the two components is reported.
##              * Direction = the response group with the higher mean
##                abundance.  Ten SGBs per direction are displayed.
##
##  R PACKAGES  mixOmics, ggplot2, patchwork
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(mixOmics)
  library(patchwork)
})
report_versions(c("mixOmics", "ggplot2", "patchwork"))


## ---- 1. Baseline SGB matrix ---------------------------------------------------

baseline <- read_stool_samples(timepoint = "Before")             # 26 patients
sgb <- read_metaphlan("sgb", samples = baseline$SampleID)         # 26 x SGBs (%), zero-sum SGBs dropped
stopifnot(all(is.finite(sgb)))
cat(nrow(sgb), "samples x", ncol(sgb), "SGBs\n")


## ---- 2. Two-component PLS-DA --------------------------------------------------
## `ncomp = 2` fixes the number of latent components (only the score plot uses
## component 2).  `scale = TRUE` = unit-variance scaling of every SGB; set it
## to FALSE to let the most abundant SGBs dominate the components.

set.seed(123)                                                     # plsda() itself is deterministic
model <- plsda(sgb, baseline$Response, ncomp = 2, scale = TRUE)

var_explained <- explained_variance(sgb, model$variates$X, ncomp = 2)   # % of SGB variance per component
cat(sprintf("Explained variance: comp1 = %.2f %%, comp2 = %.2f %%\n",
            100 * var_explained[1], 100 * var_explained[2]))

scores <- data.frame(SampleID = baseline$SampleID, Response = baseline$Response,
                     Comp1 = model$variates$X[, 1], Comp2 = model$variates$X[, 2])
save_table(scores, "Fig1H_PLS_DA_scores")


## ---- 3. VIP scores and enrichment direction -----------------------------------

vip_mat <- vip(model)                                             # SGBs x components
group_means <- sapply(response_levels, function(g) colMeans(sgb[baseline$Response == g, , drop = FALSE]))

vip_scores <- data.frame(
  SGB        = rownames(vip_mat),
  VIP_comp1  = vip_mat[, 1],
  VIP_comp2  = vip_mat[, 2],
  VIP        = apply(vip_mat, 1, max),                            # <-- max over components (as submitted)
  Mean_pCR   = group_means[rownames(vip_mat), "pCR"],
  Mean_nonpCR = group_means[rownames(vip_mat), "non-pCR"],
  row.names  = NULL
)
vip_scores$Direction <- factor(ifelse(vip_scores$Mean_pCR > vip_scores$Mean_nonpCR, "pCR", "non-pCR"),
                               levels = response_levels)
vip_scores <- vip_scores[order(-vip_scores$VIP), ]
save_table(vip_scores, "Fig1H_PLS_DA_VIP_scores")

## Ten SGBs per direction (not the twenty highest overall)  <-- change n_per_group
n_per_group <- 10
selected <- do.call(rbind, lapply(response_levels, function(g) {
  head(vip_scores[vip_scores$Direction == g, ], n_per_group)
}))
selected$Signed_VIP <- ifelse(selected$Direction == "pCR", selected$VIP, -selected$VIP)

## Display labels (plotmath): italic binomial, upright strain / SGB codes.  The
## submitted panel names SGB6001 at subspecies level; the model identifier is unchanged.
selected$Label <- taxon_plotmath(selected$SGB)
selected$Label[selected$SGB == "Fusobacterium_nucleatum|SGB6001"] <-
  'italic("Fusobacterium nucleatum")~"subsp."~italic("polymorphum")*"|SGB6001"'
selected <- selected[order(selected$Signed_VIP), ]               # most pCR-enriched at the top
selected$SGB <- factor(selected$SGB, levels = selected$SGB)
axis_labels <- setNames(selected$Label, selected$SGB)
print(selected[, c("SGB", "VIP", "Direction")], row.names = FALSE)


## ---- 4. Figure 1H -------------------------------------------------------------
## Design (as submitted): bars extend to the right for pCR-enriched and to the
## left for non-pCR-enriched SGBs; the axis shows |VIP|; the VIP value is
## printed in white inside each bar; the score plot sits in the lower-right
## corner without axis text.

x_max <- max(abs(selected$Signed_VIP)) * 1.05

p_vip <- ggplot(selected, aes(Signed_VIP, SGB, fill = Direction)) +
  geom_col(width = 0.8) +
  geom_text(aes(x = Signed_VIP * 0.85, label = sprintf("%.1f", VIP)),
            colour = "white", size = 3.2, fontface = "bold") +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  ## direction arrows above the bars
  annotate("segment", x = 0.1, xend = x_max, y = nrow(selected) + 1, yend = nrow(selected) + 1,
           colour = response_colors[["pCR"]], arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("segment", x = -0.1, xend = -x_max * 0.7, y = nrow(selected) + 1, yend = nrow(selected) + 1,
           colour = response_colors[["non-pCR"]], arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("text", x = 0.15, y = nrow(selected) + 1, label = "pCR", hjust = 0, vjust = -0.4,
           colour = response_colors[["pCR"]], fontface = "bold", size = 3.4) +
  annotate("text", x = -0.15, y = nrow(selected) + 1, label = "non-pCR", hjust = 1, vjust = -0.4,
           colour = response_colors[["non-pCR"]], fontface = "bold", size = 3.4) +
  scale_fill_manual(values = response_colors, guide = "none") +
  scale_x_continuous(labels = abs, expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_discrete(labels = function(b) parse(text = axis_labels[b]), expand = expansion(add = c(0.6, 1.6))) +
  labs(x = "VIP score", y = NULL) +
  theme_tnt(base_size = 11, legend = "none") +
  theme(axis.text.y = element_text(size = 8.5), axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

## Inset: PLS-DA score plot (component 1 vs 2), no axis text
p_scores <- ggplot(scores, aes(Comp1, Comp2, colour = Response)) +
  geom_point(size = 1.6) +
  scale_colour_manual(values = response_colors, guide = "none") +
  labs(title = "PLS-DA", x = NULL, y = NULL) +
  coord_fixed() +
  theme_tnt(base_size = 9) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.line = element_blank(), plot.title = element_text(size = 10),
        plot.background = element_rect(fill = "white", colour = NA))

p_1h <- p_vip + inset_element(p_scores, left = 0.60, bottom = 0.02, right = 0.99, top = 0.47)
save_panel(p_1h, "Fig1H_SGB_PLS_DA_VIP", width = 6.4, height = 5.2)

message("Done: 04_Fig1H_SGB_PLS_DA.R")
