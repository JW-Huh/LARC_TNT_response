## =============================================================================
##  02_Fig1D_SFig1C_Taxonomic_composition.R
##  ---------------------------------------------------------------------------
##  PANELS      Figure 1D              family-level composition, group means
##              Supplementary Fig. 1C  family-level composition, each baseline sample
##
##  PURPOSE     Show the ten most abundant bacterial families at baseline as
##              stacked relative-abundance bars, first averaged by response
##              group (1D) and then for every individual sample (S1C).
##
##  INPUT       input/microbiome/metaphlan_family.csv   MetaPhlAn 4 family
##                                                      relative abundance (%)
##              input/metadata/stool_samples.csv, input/metadata/patients.csv
##
##  OUTPUT      output/figures/Fig1D_family_composition_groups.{svg,png}
##              output/figures/SFig1C_family_composition_samples.{svg,png}
##              output/tables/Fig1D_family_group_means.csv
##
##  METHODS     Families are ranked by their summed abundance over the 26
##              baseline samples.  Figure 1D shows the top-10 families; all
##              remaining families are pooled into "Others" (= 100 - top-10 sum).
##              S1C shows the same ten families for every sample.  Bars are
##              drawn with position = "fill" so each bar sums to 1 (MetaPhlAn
##              profiles already sum to 100 %).  Samples in S1C are ordered by
##              decreasing Lachnospiraceae abundance within each group.
##              NOTE  The submitted S1C legend lists three more families
##              (Enterobacteriaceae, Lactobacillaceae, Selenomonadaceae).  They
##              belong to the after-RT top-10 and were merged into the legend
##              during figure assembly; no bar segment of S1C uses them.  Set
##              `s1c_legend_extra <- TRUE` below to add those legend entries.
##
##  R PACKAGES  ggplot2, dplyr, tidyr
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
report_versions(c("ggplot2", "dplyr", "tidyr"))


## ---- 1. Baseline family table --------------------------------------------------

baseline <- read_stool_samples(timepoint = "Before")
fam <- read_metaphlan("family", samples = baseline$SampleID)     # 26 samples x families (%)
fam_long <- data.frame(SampleID = rownames(fam), fam, check.names = FALSE) %>%
  pivot_longer(-SampleID, names_to = "Family", values_to = "Abundance") %>%
  left_join(baseline[, c("SampleID", "Response")], by = "SampleID")


## ---- 2. Top-10 families (summed over the 26 baseline samples) --------------------

n_top <- 10                                                              # <-- number of named families
top_families <- names(sort(colSums(fam), decreasing = TRUE))[1:n_top]
cat("Top-10 families (pooled baseline):", paste(top_families, collapse = ", "), "\n")
stopifnot(all(top_families %in% names(family_colors)))

## Per-group ranking, printed for reference (the per-group top-10 differ only in
## the lowest-ranked members; the panels use the pooled ranking above)
for (g in response_levels) {
  s <- sort(colSums(fam[baseline$Response == g, ]), decreasing = TRUE)[1:n_top]
  cat(g, "top-10:", paste(names(s), collapse = ", "), "\n")
}

collapse_others <- function(long, keep) {
  long %>%
    mutate(Family = ifelse(Family %in% keep, Family, "Others")) %>%
    group_by(SampleID, Response, Family) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    mutate(Family = factor(Family, levels = rev(c(keep, "Others"))))
}


## ---- 3. Figure 1D: response-group means ------------------------------------------

d1 <- collapse_others(fam_long, top_families) %>%
  group_by(Response, Family) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")
save_table(d1, "Fig1D_family_group_means")

p_1d <- ggplot(d1, aes(Response, Abundance, fill = Family)) +
  geom_col(position = "fill", width = 0.75, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = family_colors, breaks = rev(levels(d1$Family))) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance", title = "Baseline") +
  theme_tnt(base_size = 11, legend = "right") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = rel(1.1)))
save_panel(p_1d, "Fig1D_family_composition_groups", width = 3.8, height = 4.2)


## ---- 4. Supplementary Figure 1C: every baseline sample ---------------------------

s1c_legend_extra <- FALSE      # TRUE = add the three after-RT families to the legend (submitted legend)
legend_families <- c(top_families, if (s1c_legend_extra) c("Enterobacteriaceae", "Lactobacillaceae", "Selenomonadaceae"))

d2 <- collapse_others(fam_long, top_families)
d2$Family <- factor(as.character(d2$Family), levels = rev(c(legend_families, "Others")))
sample_order <- d2 %>%
  filter(Family == "Lachnospiraceae") %>%
  arrange(Response, desc(Abundance)) %>%
  pull(SampleID)
d2$SampleID <- factor(d2$SampleID, levels = sample_order)

p_s1c <- ggplot(d2, aes(SampleID, Abundance, fill = Family)) +
  geom_col(position = "fill", width = 0.9, colour = "white", linewidth = 0.2) +
  facet_grid(~Response, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = family_colors, breaks = rev(levels(d2$Family)), drop = FALSE,
                    labels = function(x) gsub("_", " ", x)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Baseline samples", y = "Relative abundance") +
  theme_tnt(base_size = 11, legend = "right") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text = element_text(size = rel(1.1), face = "plain"),
        panel.spacing.x = unit(1, "lines"))
save_panel(p_s1c, "SFig1C_family_composition_samples", width = 8, height = 4)

message("Done: 02_Fig1D_SFig1C_Taxonomic_composition.R")
