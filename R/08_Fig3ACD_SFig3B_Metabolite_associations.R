## =============================================================================
##  08_Fig3ACD_SFig3B_Metabolite_associations.R
##  ---------------------------------------------------------------------------
##  PANELS      Figure 3A   Bray-Curtis PCoA of the baseline fecal metabolome with
##                          fitted metabolite / species vectors and PERMANOVA
##              Figure 3C   Spearman correlations between 24 SCFA-related species
##                          and fecal SCFA concentrations (bubble matrix)
##              Figure 3D   E. rectale acetyl-CoA C-acetyltransferase (EC 2.3.1.9)
##                          gene abundance versus fecal butyrate
##              Suppl. Fig. 3B   twelve bile-acid ratio indices by response
##
##  PURPOSE     Link the metabolome to the metagenome of the same 20 baseline
##              stool samples: is the overall metabolite profile response-
##              associated, which SCFA producers track SCFA levels, and does a
##              species-resolved butyrate-pathway gene track butyrate?
##
##  INPUT       input/metabolome/fecal_metabolites.csv          40 metabolites (umol/g)
##              input/microbiome/metaphlan_species.csv          species abundance (%)
##              input/microbiome/humann_ec_2.3.1.9_species_stratified.tsv
##                   HUMAnN 3 EC 2.3.1.9 abundance (CPM) split by contributing species
##              input/metadata/stool_samples.csv, patients.csv
##
##  OUTPUT      output/figures/Fig3A_metabolome_PCoA, Fig3C_SCFA_species_correlations,
##              Fig3D_Erectale_EC2319_butyrate, SFig3B_bile_acid_indices  (.svg/.png)
##              output/tables/Fig3A_permanova.csv, Fig3A_envfit_vectors.csv,
##              Fig3C_SCFA_species_correlations.csv, Fig3D_EC2319_tests.csv,
##              SFig3B_bile_acid_indices.csv
##
##  METHODS
##   * Metabolome ordination (3A): metabolites detected in >= 20 % of samples
##     with non-zero variance (36), raw concentrations (no standardisation),
##     Bray-Curtis dissimilarity (vegan 2.6-6.1), classical PCoA (cmdscale,
##     k = 2); axis % = positive eigenvalue / sum of positive eigenvalues.
##     PERMANOVA: adonis2, 9,999 permutations, seed 123.  Arrows: envfit()
##     (9,999 permutations, seed 42) of seven metabolites and six species onto
##     the two PCo axes, scaled x 0.6 for display.
##   * Species-SCFA correlations (3C): Spearman rho / P (normal approximation)
##     over the 20 patients for a prespecified display list of 24 species
##     (not a new significance screen) and six SCFAs plus their sum.
##   * Enzyme (3D): E. rectale-assigned EC 2.3.1.9 CPM summed over stratified
##     rows; Spearman correlation with butyrate; Wilcoxon test (default exact
##     rule) of the enzyme abundance between response groups.
##   * Bile-acid indices (S3B): log2((sum of numerator acids + 1) /
##     (sum of denominator acids + 1)), definitions in Supplementary Table 4;
##     TUDCA/THDCA is a single combined measurement.  Wilcoxon rank-sum test.
##
##  R PACKAGES  vegan, ggplot2, ggrepel, ggbeeswarm, patchwork, dplyr
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(ggrepel)
  library(ggbeeswarm)
  library(patchwork)
})
report_versions(c("vegan", "ggplot2", "ggrepel", "ggbeeswarm", "patchwork"))


## ---- 1. Baseline metabolome matched to the same stool metagenomes ---------------------

metab <- read_metabolites(timepoint = "Before")                       # 20 patients
metabolites <- metabolite_columns()
m <- as.matrix(metab[, metabolites]); rownames(m) <- metab$SampleID

retained <- colMeans(m > 0) >= 0.20 & apply(m, 2, sd) > 0              # <-- detection / variance filter
m_ord <- m[, retained]
cat("Metabolites in the ordination:", ncol(m_ord), "of", ncol(m), "\n")

species <- read_metaphlan("species", samples = metab$SampleID, drop_empty = FALSE)   # 20 x species (%)


## ---- 2. Figure 3A: PCoA, PERMANOVA and fitted vectors ---------------------------------

d_metab <- vegdist(m_ord, method = "bray")
pcoa <- cmdscale(d_metab, k = 2, eig = TRUE)
axis_pct <- 100 * pmax(pcoa$eig, 0) / sum(pmax(pcoa$eig, 0))

set.seed(123)                                                          # <-- permutation seed
permanova <- adonis2(d_metab ~ Response, data = metab, permutations = 9999)
print(permanova)
save_table(data.frame(term = rownames(permanova), as.data.frame(permanova)), "Fig3A_permanova")

## Features projected as arrows  <-- edit these two vectors to show other features
vector_metabolites <- c("Butyrate", "Acetate", "Propionate", "Indole_acetic_acid",
                        "Indole_lactic_acid", "Indolepropionic_acid", "Nicotinic_acid")
vector_species <- c("Coprococcus_comes", "Eubacterium_rectale", "Dorea_longicatena",
                    "Phocaeicola_vulgatus", "Roseburia_inulinivorans", "Eggerthella_lenta")
stopifnot(all(vector_species %in% colnames(species)))

set.seed(42)
fit <- envfit(pcoa$points, cbind(m[, vector_metabolites], species[, vector_species]),
              permutations = 9999, w = rep(1, nrow(m)))
vectors <- as.data.frame(scores(fit, display = "vectors"))
names(vectors) <- c("PCo1", "PCo2")
vectors$Feature <- rownames(vectors)
vectors$r2 <- fit$vectors$r; vectors$P <- fit$vectors$pvals
vectors$Group <- case_when(vectors$Feature %in% c("Butyrate", "Acetate", "Propionate") ~ "SCFA",
                           vectors$Feature %in% vector_species ~ "Species", TRUE ~ "Tryptophan metabolites")
vectors$Label <- gsub("_", " ", vectors$Feature)
vectors$Label[vectors$Feature == "Indolepropionic_acid"] <- "Indole propionic acid"
save_table(vectors, "Fig3A_envfit_vectors")

arrow_scale <- 0.6                                                     # <-- arrow length relative to unit scores
vectors <- vectors %>% mutate(x_end = PCo1 * arrow_scale, y_end = PCo2 * arrow_scale)

scores_df <- data.frame(pcoa$points, Response = metab$Response); names(scores_df)[1:2] <- c("PCo1", "PCo2")
label_colors <- c("SCFA" = "#6EA87D", "Tryptophan metabolites" = "#7181A5", "Species" = "#C28968")
permanova_label <- sprintf('paste("PERMANOVA: ", R^2, " = %.3f, ", italic(p), " = %.3f")',
                           permanova$R2[1], permanova$`Pr(>F)`[1])

p_ord <- ggplot(scores_df, aes(PCo1, PCo2)) +
  geom_segment(data = vectors, aes(x = 0, y = 0, xend = x_end, yend = y_end), colour = "grey45",
               arrow = arrow(length = unit(0.16, "cm")), linewidth = 0.4) +
  geom_point(aes(colour = Response), size = 2.4) +
  geom_text_repel(data = vectors, aes(x_end, y_end, label = Label, colour = Group,
                                      fontface = ifelse(Group == "Species", "italic", "plain")),
                  size = 3, seed = 42, max.overlaps = Inf, segment.colour = NA, show.legend = FALSE,
                  box.padding = 0.25, point.padding = 0.2) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8, size = 3.4, label = permanova_label, parse = TRUE) +
  scale_colour_manual(values = c(response_colors, label_colors), guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  scale_y_continuous(expand = expansion(mult = 0.12)) +
  labs(x = sprintf("PCo1 (%.1f%%)", axis_pct[1]), y = sprintf("PCo2 (%.1f%%)", axis_pct[2])) +
  theme_tnt(base_size = 11, legend = "none") +
  theme(aspect.ratio = 1,                                               # square panel as in the submitted figure
        axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6))

## hand-made legend strip: coloured dots for the groups, coloured "a" for the label colours
legend_df <- data.frame(x = c(0, 0.9, 2.2, 3.4), y = 0,
                        key = c("pCR", "non-pCR", "SCFA", "Tryptophan metabolites"),
                        type = c("point", "point", "text", "text"))
p_key <- ggplot(legend_df, aes(x, y)) +
  geom_point(data = filter(legend_df, type == "point"), aes(colour = key), size = 3) +
  geom_text(data = filter(legend_df, type == "text"), aes(colour = key), label = "a", size = 3.6) +
  geom_text(aes(x = x + 0.12, label = key), hjust = 0, size = 3.4) +
  scale_colour_manual(values = c(response_colors, label_colors), guide = "none") +
  coord_cartesian(xlim = c(-0.2, 5.6), ylim = c(-0.5, 0.5)) +
  theme_void()

p_3a <- p_key / p_ord + plot_layout(heights = c(0.06, 1))
save_panel(p_3a, "Fig3A_metabolome_PCoA", width = 4.6, height = 4.8)


## ---- 3. Figure 3C: SCFA-related species versus SCFA concentrations --------------------
## A prespecified display list of 24 species (butyrate / propionate producers
## and cross-feeders discussed in the text), not a new significance screen.

scfa_species <- c(
  "Eubacterium_rectale", "Coprococcus_comes", "Megasphaera_elsdenii", "Ruminococcus_callidus",
  "Ruminococcus_bromii", "Coprococcus_catus", "Anaerostipes_caccae", "Roseburia_inulinivorans",
  "Roseburia_intestinalis", "Phascolarctobacterium_succinatutens", "Coprococcus_eutactus",
  "Faecalibacterium_prausnitzii", "Blautia_hydrogenotrophica", "Anaerostipes_butyraticus",
  "Anaerobutyricum_hallii", "Eubacterium_ramulus", "Blautia_obeum", "Anaerobutyricum_soehngenii",
  "Roseburia_hominis", "Eubacterium_limosum", "Roseburia_faecis", "Anaerostipes_hadrus",
  "Akkermansia_muciniphila", "Anaerotruncus_colihominis"
)
scfas <- c("Isovalerate", "Isobutyrate", "Valerate", "Propionate", "Acetate", "Butyrate")
stopifnot(all(scfa_species %in% colnames(species)))
metab$Total_SCFA <- rowSums(metab[, scfas])

cors <- expand.grid(Species = scfa_species, Metabolite = c(scfas, "Total_SCFA"), stringsAsFactors = FALSE)
cors[, c("Rho", "P")] <- t(mapply(function(sp, met) {
  ct <- suppressWarnings(cor.test(species[, sp], metab[[met]], method = "spearman", exact = FALSE))
  c(unname(ct$estimate), ct$p.value)
}, cors$Species, cors$Metabolite))
save_table(cors, "Fig3C_SCFA_species_correlations")

cors$Species <- factor(cors$Species, levels = scfa_species)
cors$Metabolite <- factor(cors$Metabolite, levels = c(scfas, "Total_SCFA"))
## two-line axis labels (genus / species) as in the submitted panel; all 24 names
## are Latin binomials, so the whole label is set in italics
species_labels <- setNames(sub("_", "\n", scfa_species), scfa_species)

p_3c <- ggplot(cors, aes(Species, Metabolite)) +
  geom_point(aes(fill = Rho, size = -log10(P)), shape = 21, colour = "black", stroke = 0.4) +
  scale_fill_gradient2(low = "#D1606E", mid = "white", high = "#4C7FB6", midpoint = 0,
                       limits = c(-0.5, 0.6), oob = scales::squish, breaks = c(-0.25, 0, 0.25, 0.5),
                       name = expression("Spearman's " * rho)) +
  scale_size_continuous(range = c(1.2, 8), breaks = c(0, 0.5, 1, 1.5, 2), limits = c(0, 2.3),
                        name = expression(-Log[10](italic(p)))) +
  scale_x_discrete(labels = function(b) species_labels[b]) +
  scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
  guides(fill = guide_colourbar(order = 1, barwidth = unit(3, "cm"), barheight = unit(0.3, "cm")),
         size = guide_legend(order = 2, nrow = 1)) +
  labs(x = NULL, y = NULL) +
  theme_tnt(base_size = 9, legend = "top") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8.5, face = "italic", lineheight = 0.9),
        axis.line = element_blank(),
        axis.text.y = element_text(size = 9),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        legend.box = "horizontal", legend.title = element_text(size = 9))
save_panel(p_3c, "Fig3C_SCFA_species_correlations", width = 8.4, height = 4.4)


## ---- 4. Figure 3D: E. rectale EC 2.3.1.9 abundance versus fecal butyrate ---------------

ec <- read.delim(input_file("microbiome", "humann_ec_2.3.1.9_species_stratified.tsv"), check.names = FALSE, quote = "")
rect_rows <- grepl("s__Eubacterium_rectale$", ec$Taxon)
stopifnot(sum(rect_rows) >= 1)
metab$EC_2319_Erectale <- colSums(as.matrix(ec[rect_rows, metab$SampleID, drop = FALSE]))   # CPM

rho_test <- cor.test(metab$EC_2319_Erectale, metab$Butyrate, method = "spearman", exact = FALSE)
group_p  <- wilcox_p(metab$EC_2319_Erectale, metab$Response)
cat(sprintf("EC 2.3.1.9 (E. rectale) vs butyrate: rho = %.3f, P = %.4f; group Wilcoxon P = %.4f\n",
            rho_test$estimate, rho_test$p.value, group_p))
save_table(data.frame(Test = c("Spearman rho", "Spearman P", "Wilcoxon P (pCR vs non-pCR)"),
                      Value = c(rho_test$estimate, rho_test$p.value, group_p)), "Fig3D_EC2319_tests")

p_scatter <- ggplot(metab, aes(EC_2319_Erectale, Butyrate)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2E9BDA", fill = "#BEE5F2", alpha = 0.6, linewidth = 0.6) +
  geom_point(aes(colour = Response), size = 2.6) +
  annotate("text", x = mean(range(metab$EC_2319_Erectale)), y = Inf, vjust = 1.8, size = 3.4, parse = TRUE,
           label = sprintf('paste("Spearman\'s ", rho, " = %.2f , ", italic(p), " = %.2f")', rho_test$estimate, rho_test$p.value)) +
  scale_colour_manual(values = response_colors, guide = "none") +
  labs(x = "Acetyl-CoA C-acetyltransferase",
       caption = expression("(" * italic("E. rectale") * ")"),     # second line of the x title (plotmath cannot break lines)
       y = expression("Fecal butyrate (" * mu * "mol/g)")) +
  theme_tnt(base_size = 10, legend = "none") +
  theme(plot.caption = element_text(hjust = 0.5, size = 10, colour = "black", margin = margin(t = 1)),
        axis.title.x = element_text(margin = margin(t = 4, b = 0)))

## inset: enzyme abundance by response group (horizontal boxes)
p_box <- ggplot(metab, aes(EC_2319_Erectale, Response, colour = Response, fill = Response)) +
  geom_boxplot(width = 0.5, alpha = 0.35, outlier.shape = NA, linewidth = 0.5) +
  annotate("segment", x = Inf, xend = Inf, y = 1, yend = 2, linewidth = 0.4) +
  annotate("text", x = Inf, y = 1.5, label = formatC(group_p, format = "f", digits = 3), angle = 90, vjust = -0.4, size = 2.6) +
  scale_colour_manual(values = response_colors, guide = "none") +
  scale_fill_manual(values = response_colors, guide = "none") +
  scale_y_discrete(limits = response_levels) +                          # pCR bottom, non-pCR top (as submitted)
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_tnt(base_size = 8, legend = "none") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        plot.background = element_blank(), panel.background = element_blank(), plot.margin = margin(2, 12, 2, 2))

p_3d <- p_scatter + inset_element(p_box, left = 0.45, bottom = 0.05, right = 0.98, top = 0.28)
save_panel(p_3d, "Fig3D_Erectale_EC2319_butyrate", width = 4.4, height = 4.2)


## ---- 5. Supplementary Figure 3B: bile-acid ratio indices --------------------------------

ba_columns <- c(
  CA = "Cholic_acid", CDCA = "Chenodeoxycholic_acid", GCA = "Glycocholic_acid",
  GCDCA = "Glycochenodeoxycholic_acid", TCA = "Taurocholic_acid", TCDCA = "Taurochenodeoxycholic_acid",
  DCA = "Deoxycholic_acid", HDCA = "Hyodeoxycholic_acid", LCA = "Lithocholic_acid",
  UDCA = "Ursodeoxycholic_acid", GDCA = "Glycodeoxycholic_acid", GLCA = "Glycolithocholic_acid",
  GUDCA = "Glycoursodeoxycholic_acid", TDCA = "Taurodeoxycholic_acid", TLCA = "Taurolithocholic_acid",
  TUDCA_THDCA = "Tauroursodeoxycholic_acid_Taurohyodeoxycholic_acid"
)
ba <- as.data.frame(m[, ba_columns]); names(ba) <- names(ba_columns)

## Index definitions (Supplementary Table 4): numerator / denominator acid sets
bile_indices <- list(
  "Secondary/Primary"        = list(num = c("DCA", "LCA", "UDCA", "HDCA", "GDCA", "TDCA", "GLCA", "TLCA", "GUDCA", "TUDCA_THDCA"),
                                    den = c("CA", "CDCA", "GCA", "TCA", "GCDCA", "TCDCA")),
  "Conjugated/Unconjugated"  = list(num = c("GCA", "TCA", "GCDCA", "TCDCA", "GDCA", "TDCA", "GLCA", "TLCA", "GUDCA", "TUDCA_THDCA"),
                                    den = c("CA", "CDCA", "DCA", "LCA", "UDCA", "HDCA")),
  "Glycine/Taurine"          = list(num = c("GCA", "GCDCA", "GDCA", "GLCA", "GUDCA"),
                                    den = c("TCA", "TCDCA", "TDCA", "TLCA", "TUDCA_THDCA")),
  "CA-derived/CDCA-derived"  = list(num = c("CA", "DCA", "GCA", "TCA", "GDCA", "TDCA"),
                                    den = c("CDCA", "LCA", "UDCA", "HDCA", "GCDCA", "TCDCA", "GLCA", "TLCA", "GUDCA", "TUDCA_THDCA")),
  "Hydrophobic/Hydrophilic"  = list(num = c("CDCA", "DCA", "LCA", "GCDCA", "TCDCA", "GDCA", "TDCA", "GLCA", "TLCA"),
                                    den = c("CA", "UDCA", "GCA", "TCA", "GUDCA", "TUDCA_THDCA")),
  "Cytotoxic/Protective"     = list(num = c("DCA", "LCA"), den = "UDCA"),
  "Total 7α"            = list(num = c("DCA", "LCA"), den = c("CA", "CDCA")),
  "CA-derived 7α"       = list(num = "DCA", den = "CA"),
  "CDCA-derived 7α"     = list(num = "LCA", den = "CDCA"),
  "Total BSH"                = list(num = c("CA", "CDCA"), den = c("GCA", "TCA", "GCDCA", "TCDCA")),
  "CA-derived BSH"           = list(num = "CA", den = c("GCA", "TCA")),
  "CDCA-derived BSH"         = list(num = "CDCA", den = c("GCDCA", "TCDCA"))
)

index_values <- sapply(bile_indices, function(idx) {
  log2((rowSums(ba[, idx$num, drop = FALSE]) + 1) / (rowSums(ba[, idx$den, drop = FALSE]) + 1))
})
index_tests <- data.frame(Index = names(bile_indices),
                          P = apply(index_values, 2, function(v) wilcox_p(v, metab$Response)))
print(index_tests, row.names = FALSE, digits = 3)
save_table(cbind(metab[, c("SubjectID", "SampleID", "Response")], index_values), "SFig3B_bile_acid_indices")

## Design (as submitted): jittered points, black median bar, P above, blue
## median-difference bracket (as in Figure 1B)
index_panel <- function(name, show_x = TRUE) {
  d <- data.frame(Response = metab$Response, Value = index_values[, name])
  p <- index_tests$P[index_tests$Index == name]
  med <- tapply(d$Value, d$Response, median)
  y_top <- max(d$Value) + 0.08 * diff(range(d$Value))
  ggplot(d, aes(Response, Value)) +
    geom_quasirandom(aes(colour = Response), width = 0.2, size = 1.6) +
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.4, linewidth = 0.4, colour = "black") +                 # median bar
    annotate("segment", x = 1, xend = 2, y = y_top, yend = y_top, linewidth = 0.35) +
    annotate("text", x = 1.5, y = y_top, vjust = -0.5, size = 2.7, label = formatC(p, format = "g", digits = 2)) +
    annotate("segment", x = 1.35, xend = 1.45, y = med["pCR"], yend = med["pCR"], colour = "#2988CA", linewidth = 0.35) +
    annotate("segment", x = 1.45, xend = 1.55, y = med["non-pCR"], yend = med["non-pCR"], colour = "#2988CA", linewidth = 0.35) +
    annotate("segment", x = 1.45, xend = 1.45, y = med["pCR"], yend = med["non-pCR"], colour = "#2988CA", linewidth = 0.35) +
    annotate("text", x = 1.45, y = mean(med), label = as.character(round(med["pCR"] - med["non-pCR"], 2)),   # median difference pCR - non-pCR
             colour = "#2988CA", size = 2.3, hjust = -0.2, vjust = -0.5) +
    scale_colour_manual(values = response_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.2))) +
    labs(x = NULL, y = NULL, title = name) +
    theme_tnt(base_size = 8, legend = "none") +
    theme(plot.title = element_text(size = 8),
          axis.text.x = if (show_x) element_text(size = 8) else element_blank())
}
## x tick labels only on the bottom row, as in the submitted figure
p_s3b <- wrap_plots(lapply(seq_along(bile_indices), function(i) index_panel(names(bile_indices)[i], show_x = i > length(bile_indices) - 3)), ncol = 3) +
  plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
p_s3b <- wrap_elements(p_s3b) + labs(tag = expression(Log[2] * "(Bile Acids Ratio)")) +
  theme(plot.tag = element_text(angle = 90, size = 10), plot.tag.position = "left")
save_panel(p_s3b, "SFig3B_bile_acid_indices", width = 5.6, height = 8.4)

message("Done: 08_Fig3ACD_SFig3B_Metabolite_associations.R")
