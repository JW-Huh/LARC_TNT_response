## =============================================================================
##  01_Fig1BC_SFig1ABD_Diversity.R
##  ---------------------------------------------------------------------------
##  PANELS      Figure 1B  Pielou's evenness (baseline, pCR vs non-pCR)
##              Figure 1C  Inverse Simpson diversity
##              Supplementary Figure 1A  Shannon index
##              Supplementary Figure 1B  Observed SGB richness
##              Supplementary Figure 1D  Bray-Curtis PCoA + PERMANOVA
##
##  PURPOSE     Compare baseline gut-microbiome alpha diversity between
##              pathological complete responders (pCR, TRG 0) and non-responders
##              (non-pCR, TRG 1-3), and test whether overall community
##              composition differs between the groups.
##
##  INPUT       input/microbiome/metaphlan_sgb.csv   MetaPhlAn 4 SGB relative
##                                                   abundance (%), 42 samples
##              input/metadata/stool_samples.csv, input/metadata/patients.csv
##
##  OUTPUT      output/figures/Fig1B_Evenness, Fig1C_Inverse_Simpson,
##              SFig1A_Shannon, SFig1B_Richness, SFig1D_PCoA_Bray_Curtis  (.svg/.png)
##              output/tables/Fig1BC_SFig1AB_alpha_diversity.csv
##              output/tables/SFig1D_permanova.csv
##
##  METHODS     * Alpha diversity (vegan 2.6-6.1): richness = number of SGBs with
##                abundance > 0; Shannon H = -sum(p_i log p_i); inverse Simpson
##                = 1 / sum(p_i^2); Pielou's evenness J = H / ln(richness).
##                Computed on relative abundances (%) without rarefaction.
##              * Group comparison: two-sided Wilcoxon rank-sum test (R default
##                settings: exact P when there are no ties, n < 50).
##              * Beta diversity: Bray-Curtis dissimilarity of SGB profiles,
##                principal-coordinate analysis (cmdscale, 10 axes; the axis %
##                is eigenvalue / sum of all eigenvalues as in the submission),
##                PERMANOVA with vegan::adonis2 (999 permutations, seed 123).
##
##  R PACKAGES  vegan, ggplot2, ggbeeswarm
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(vegan)
  library(ggbeeswarm)
})
report_versions(c("vegan", "ggplot2", "ggbeeswarm"))


## ---- 1. Baseline SGB profiles -------------------------------------------------

baseline <- read_stool_samples(timepoint = "Before")            # 26 patients
sgb <- read_metaphlan("sgb", samples = baseline$SampleID)        # 26 x SGBs (% abundance)
stopifnot(nrow(sgb) == 26, all(rowSums(sgb) > 0))


## ---- 2. Alpha diversity ------------------------------------------------------

alpha <- data.frame(
  SampleID   = baseline$SampleID,
  SubjectID  = baseline$SubjectID,
  Response   = baseline$Response,
  Richness   = specnumber(sgb),                    # observed SGBs
  Shannon    = diversity(sgb, index = "shannon"),
  InvSimpson = diversity(sgb, index = "invsimpson")
)
alpha$Evenness <- alpha$Shannon / log(alpha$Richness)   # Pielou's J

## Wilcoxon rank-sum tests (two-sided).  `exact = NULL` = R's default rule.
alpha_tests <- data.frame(
  Index = c("Evenness", "InvSimpson", "Shannon", "Richness"),
  P = sapply(c("Evenness", "InvSimpson", "Shannon", "Richness"),
             function(v) wilcox_p(alpha[[v]], alpha$Response))
)
print(alpha_tests)
save_table(alpha, "Fig1BC_SFig1AB_alpha_diversity")
save_table(alpha_tests, "Fig1BC_SFig1AB_alpha_diversity_tests")


## ---- 3. Violin + box + beeswarm panels ----------------------------------------
## Design (as submitted): translucent violin, whisker-less box, individual
## points, a P-value bracket on top and a blue bracket showing the difference
## in group medians (#2988CA).  Change `panel_size` to resize every panel.

panel_size <- c(width = 2.3, height = 4.2)

alpha_panel <- function(index, title, digits_delta = 3) {
  d <- data.frame(Response = alpha$Response, Value = alpha[[index]])
  p <- alpha_tests$P[alpha_tests$Index == index]

  med <- tapply(d$Value, d$Response, median)
  delta <- unname(abs(med["pCR"] - med["non-pCR"]))
  y_range <- range(d$Value)
  y_top <- y_range[2] + 0.06 * diff(y_range)        # P-value bracket height
  y_show <- y_range + c(-0.10, 0.22) * diff(y_range) # visible range (violin tails are clipped)

  ggplot(d, aes(Response, Value)) +
    geom_violin(aes(fill = Response), colour = NA, alpha = 0.3, width = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, coef = 0, linewidth = 0.5,
                 fill = NA, colour = "grey30") +
    geom_beeswarm(aes(colour = Response), size = 2.3, cex = 2.5) +
    ## P-value bracket
    annotate("segment", x = 1, xend = 2, y = y_top, yend = y_top, linewidth = 0.4) +
    annotate("text", x = 1.5, y = y_top, vjust = -0.6, size = 3.6,
             label = paste0("italic(p) == ", formatC(p, format = "f", digits = 3)), parse = TRUE) +
    ## median-difference bracket (blue)
    annotate("segment", x = 1.35, xend = 1.45, y = med["pCR"], yend = med["pCR"], colour = "#2988CA", linewidth = 0.4) +
    annotate("segment", x = 1.45, xend = 1.55, y = med["non-pCR"], yend = med["non-pCR"], colour = "#2988CA", linewidth = 0.4) +
    annotate("segment", x = 1.45, xend = 1.45, y = med["pCR"], yend = med["non-pCR"], colour = "#2988CA", linewidth = 0.4) +
    annotate("text", x = 1.45, y = mean(med), label = format(round(delta, digits_delta), nsmall = 0),
             colour = "#2988CA", size = 3, hjust = -0.15, vjust = -0.6) +
    scale_fill_manual(values = response_colors) +
    scale_colour_manual(values = response_colors) +
    coord_cartesian(ylim = y_show) +
    labs(x = NULL, y = NULL, title = title) +
    theme_tnt(base_size = 11, legend = "none") +
    theme(axis.text.x = element_text(size = rel(1.1)))
}

save_panel(alpha_panel("Evenness",   "Pielou's Evenness", 3), "Fig1B_Evenness",         panel_size[1], panel_size[2])
save_panel(alpha_panel("InvSimpson", "Inverse Simpson",   2), "Fig1C_Inverse_Simpson",  panel_size[1], panel_size[2])
save_panel(alpha_panel("Shannon",    "Shannon Index",     2), "SFig1A_Shannon",         panel_size[1], panel_size[2])
save_panel(alpha_panel("Richness",   "Observed SGBs",     0), "SFig1B_Richness",        panel_size[1], panel_size[2])


## ---- 4. Supplementary Figure 1D: Bray-Curtis PCoA and PERMANOVA ---------------

bray <- vegdist(sgb, method = "bray")
pcoa <- cmdscale(bray, k = 10, eig = TRUE)

## Percent variance of each axis (denominator = sum of ALL eigenvalues, as in
## the submitted panel: 13.96 % and 10.45 %)
axis_pct <- round(100 * pcoa$eig / sum(pcoa$eig), 2)

set.seed(123)                                       # <-- permutation seed
permanova <- adonis2(bray ~ Response, data = baseline, permutations = 999)
print(permanova)
save_table(data.frame(term = rownames(permanova), as.data.frame(permanova)), "SFig1D_permanova")

scores <- data.frame(pcoa$points[, 1:2], Response = baseline$Response)
names(scores)[1:2] <- c("PCo1", "PCo2")
hulls <- do.call(rbind, lapply(split(scores, scores$Response), function(d) d[chull(d$PCo1, d$PCo2), ]))

p_pcoa <- ggplot(scores, aes(PCo1, PCo2)) +
  geom_polygon(data = hulls, aes(fill = Response, colour = Response), alpha = 0.08, linewidth = 0.8) +
  geom_point(aes(fill = Response), shape = 21, size = 3.5, colour = "white", stroke = 0.3) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 3.6,
           label = sprintf("R² = %.2f %%, PERMANOVA p = %.3f",
                           100 * permanova$R2[1], permanova$`Pr(>F)`[1])) +
  scale_fill_manual(values = response_colors, name = "TRG") +
  scale_colour_manual(values = response_colors, name = "TRG") +
  coord_fixed() +
  labs(x = sprintf("PCo1 (%.2f%%)", axis_pct[1]), y = sprintf("PCo2 (%.2f%%)", axis_pct[2])) +
  theme_tnt(base_size = 11, legend = "inside") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.3))
save_panel(p_pcoa, "SFig1D_PCoA_Bray_Curtis", width = 4.2, height = 4.2)

message("Done: 01_Fig1BC_SFig1ABD_Diversity.R")
