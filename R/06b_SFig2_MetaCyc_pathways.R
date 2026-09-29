## =============================================================================
##  06b_SFig2_MetaCyc_pathways.R
##  ---------------------------------------------------------------------------
##  PANELS      Supplementary Figure 2A    heatmap of response-associated MetaCyc pathways
##              Supplementary Figure 2B/C  sucrose and rhamnose "pathway balance"
##                                         groups versus TRG score (scatter + tiles)
##              Supplementary Figure 2D-G  genus contributions to four pathways
##              Supplementary Figure 2H-M  species contributions (Bacteroides, Blautia)
##
##  PURPOSE     Describe which MetaCyc pathways differ at baseline, whether a
##              carbohydrate-pathway "balance" (biosynthesis vs degradation of
##              sucrose / rhamnose) tracks tumour regression, and which taxa
##              encode those pathways.
##
##  INPUT       input/microbiome/humann_pathway_abundance.tsv   HUMAnN 3 MetaCyc
##                    pathway abundance (unstratified and taxon-stratified rows)
##              input/microbiome/humann_pathway_coverage.tsv    matching coverage
##              input/microbiome/metaphlan_species.csv          species abundance (%)
##              input/metadata/stool_samples.csv, patients.csv
##
##  OUTPUT      output/figures/SFig2A_pathway_heatmap, SFig2B_sucrose_balance,
##              SFig2C_rhamnose_balance, SFig2D-G_<pathway>_genus_contributions,
##              SFig2H-M_<genus>_<pathway>  (.svg/.png)
##              output/tables/SFig2A_pathway_tests.csv, SFig2BC_pathway_balance.csv
##
##  METHODS
##   * Coverage mask: a pathway's abundance in a sample is set to 0 when its
##     HUMAnN coverage is < 0.5 (unlikely to be complete in that sample).
##   * Pathway screen (S2A): the figure-producing code kept pathways with
##     > 0 abundance in at least 0.2 x 26 = 5.2 of the 42 metagenomes and a mean
##     coverage >= 1/26 over the 42 metagenomes (so this is not a strictly
##     baseline-only filter - kept as run; see REVIEW notes).  Two-sided
##     Wilcoxon rank-sum test (normal approximation) on the 26 baseline samples;
##     pathways with P < 0.2 are displayed (11); significance stars in the
##     P-value annotation: ** P < 0.05, * P < 0.10 (thresholds of the submitted
##     panel and its legend).  Heatmap: log10(x + pseudocount),
##     pseudocount = half of the smallest positive value; row z-scores; rows
##     clustered with Manhattan distance / average linkage; samples clustered
##     within response group (Euclidean / complete).  "LFC" is log10 of the
##     ratio of group means (the submitted legend calls it Log2FC; the values are
##     log10 - the label here is corrected).
##   * Balance groups (S2B/C): sucrose-high = PWY-7238 (sucrose biosynthesis II)
##     > 30,000 AND degradation (PWY-621 + PWY-5384) == 0; rhamnose-high =
##     DTDPRHAMSYN-PWY >= 35,000 OR RHAMCAT-PWY <= 5,000.  TRG score compared
##     between groups by Wilcoxon rank-sum test with continuity correction.
##   * Contributions (S2D-M): taxon-stratified pathway abundance (coverage-
##     masked), summed per response group and divided by group size (group mean
##     including zeros).  Genus panels show the six largest genera; species
##     panels show every species of the genus that contributes to the pathways.
##     Three Bacteroides -> Phocaeicola renamings (vulgatus, dorei, plebeius) are
##     harmonised between HUMAnN and MetaPhlAn names.
##
##  R PACKAGES  dplyr, tidyr, ggplot2, patchwork, ComplexHeatmap, circlize
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})
report_versions(c("ggplot2", "dplyr", "tidyr", "patchwork", "ComplexHeatmap", "circlize"))
set_heatmap_options()                                   # legend font sizes (R/_common.R)


## ---- 1. HUMAnN pathway abundance and coverage ------------------------------------

baseline <- read_stool_samples(timepoint = "Before")
baseline <- baseline[order(baseline$SampleID), ]                 # lexical order, as in the original run

abundance <- read_feature_table("humann_pathway_abundance.tsv", "Pathway")     # all 42 samples
coverage  <- read_feature_table("humann_pathway_coverage.tsv",  "Pathway")
stopifnot(identical(dimnames(abundance), dimnames(coverage)), all(is.finite(abundance)))

coverage_cut <- 0.5                                              # <-- HUMAnN coverage mask
unstratified <- !grepl("|", rownames(abundance), fixed = TRUE)
raw <- abundance[unstratified, ]
cov <- coverage[unstratified, ]
base_values <- raw[, baseline$SampleID]
base_values[cov[, baseline$SampleID] < coverage_cut] <- 0

## Pathway screen (kept exactly as in the figure-producing code, see header)
retained <- rowSums(raw > 0) >= 0.2 * nrow(baseline) & rowMeans(cov) >= 1 / nrow(baseline) &
            !rownames(raw) %in% c("UNMAPPED", "UNINTEGRATED")

tests <- data.frame(
  Pathway = rownames(raw)[retained],
  P = apply(base_values[retained, ], 1, function(x) wilcox_p(x, baseline$Response, exact = FALSE))
)
display_p_cut <- 0.2                                             # <-- heatmap inclusion threshold
selected <- sort(tests$Pathway[!is.na(tests$P) & tests$P < display_p_cut])
cat("Screened pathways:", sum(retained), "; displayed (P <", display_p_cut, "):", length(selected), "\n")
save_table(tests[order(tests$P), ], "SFig2A_pathway_tests")


## ---- 2. Supplementary Figure 2A: pathway heatmap ----------------------------------

pseudo  <- min(base_values[base_values > 0]) / 2
values  <- base_values[selected, ]
z       <- t(scale(t(log10(values + pseudo))))
p_sel   <- tests$P[match(selected, tests$Pathway)]
lfc10   <- log10((rowMeans(values[, baseline$Response == "pCR"]) + pseudo) /
                 (rowMeans(values[, baseline$Response == "non-pCR"]) + pseudo))
mean_cov <- rowMeans(cov[selected, baseline$SampleID])

col_order <- order_samples_within_response(z, baseline)
meta <- baseline[match(col_order, baseline$SampleID), ]

star_cuts <- c(0.05, 0.10); star_symbols <- c("**", "*")           # <-- star thresholds (** P < 0.05, * P < 0.10)
right <- rowAnnotation(
  `p-value` = anno_simple(-log10(p_sel), col = colorRamp2(c(0.2, 1.8), c("white", "#2166AC")),
                          pch = p_stars(p_sel, star_cuts, star_symbols), pt_gp = gpar(col = "white", fontsize = 8),
                          gp = gpar(col = "grey40", lwd = 0.5), width = unit(3.5, "mm")),
  Coverage = anno_simple(mean_cov, col = colorRamp2(c(0, 1), c("white", "#B24A5A")),
                         gp = gpar(col = "grey40", lwd = 0.5), width = unit(3.5, "mm")),
  LFC = anno_simple(lfc10, col = colorRamp2(c(-2, 0, 2), c("#B2182B", "white", "#4FAE9A")),
                    gp = gpar(col = "grey40", lwd = 0.5), width = unit(3.5, "mm")),
  annotation_name_gp = gpar(fontsize = 8), annotation_name_rot = 90
)
ht <- Heatmap(
  z[, col_order], name = "z-score", col = colorRamp2(c(-2, 0, 2), zscore_colors),
  cluster_columns = FALSE, column_split = meta$Response, column_title = NULL, column_gap = unit(1, "mm"),
  clustering_distance_rows = "manhattan", clustering_method_rows = "average",
  row_dend_side = "left", row_dend_width = unit(1, "cm"),
  show_column_names = FALSE, row_names_side = "right", row_names_gp = gpar(fontsize = 8),
  row_names_max_width = unit(12, "cm"), rect_gp = gpar(col = "grey40", lwd = 0.5),
  top_annotation = response_top_annotation(meta), right_annotation = right,
  heatmap_legend_param = list(direction = "horizontal", at = -2:2, title_position = "topcenter")
)
legends <- list(
  Legend(title = expression(-Log[10](italic(p))), col_fun = colorRamp2(c(0.2, 1.8), c("white", "#2166AC")),
         at = c(0.2, 0.6, 1, 1.4, 1.8), direction = "horizontal", title_position = "topcenter"),
  Legend(title = "Coverage", col_fun = colorRamp2(c(0, 1), c("white", "#B24A5A")), at = seq(0, 1, 0.2),
         direction = "horizontal", title_position = "topcenter"),
  Legend(title = expression(Log[10]*FC), col_fun = colorRamp2(c(-2, 0, 2), c("#B2182B", "white", "#4FAE9A")),
         at = -2:2, direction = "horizontal", title_position = "topcenter")
)
dev <- open_panel("SFig2A_pathway_heatmap", width = 11, height = 4.2)
draw_panel(dev, draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                     merge_legends = TRUE, annotation_legend_list = legends))


## ---- 3. Supplementary Figure 2B/C: carbohydrate pathway balance --------------------

pathway_ids <- sub(":.*$", "", rownames(base_values))
pw <- function(id) base_values[match(id, pathway_ids), ]

balance <- data.frame(
  baseline[, c("SampleID", "SubjectID", "Response", "TRG_score")],
  Sucrose_synthesis    = pw("PWY-7238"),
  Sucrose_degradation  = pw("PWY-621") + pw("PWY-5384"),
  Rhamnose_synthesis   = pw("DTDPRHAMSYN-PWY"),
  Rhamnose_degradation = pw("RHAMCAT-PWY")
)

## Group definitions  <-- change the cut-offs here (abundance units: HUMAnN CPM-like)
balance$Sucrose_group  <- ifelse(balance$Sucrose_synthesis > 30000 & balance$Sucrose_degradation == 0, "High", "Low")
balance$Rhamnose_group <- ifelse(balance$Rhamnose_synthesis >= 35000 | balance$Rhamnose_degradation <= 5000, "High", "Low")
save_table(balance, "SFig2BC_pathway_balance")

## x axis = biosynthesis pathway, y axis = degradation pathway (as submitted)
balance_panel <- function(carb, x_lab, y_lab, group_lab, synthesis_cut, degradation_cut, name) {
  d <- balance
  d$Synthesis <- d[[paste0(carb, "_synthesis")]]; d$Degradation <- d[[paste0(carb, "_degradation")]]
  d$Group <- factor(d[[paste0(carb, "_group")]], levels = c("High", "Low"))
  p <- wilcox.test(TRG_score ~ Group, data = d, exact = FALSE, correct = TRUE)$p.value
  cat(carb, "balance: Wilcoxon P =", signif(p, 3), "\n"); print(table(d$Group, d$TRG_score))

  scatter <- ggplot(d, aes(Synthesis, Degradation, fill = Response)) +
    geom_vline(xintercept = synthesis_cut, linetype = 2, colour = "grey55") +
    geom_hline(yintercept = degradation_cut, linetype = 2, colour = "grey55") +
    geom_point(shape = 21, size = 3, colour = "white") +
    scale_fill_manual(values = response_colors, name = "TRG") +
    labs(x = x_lab, y = y_lab) +
    theme_tnt(base_size = 9, legend = "none")

  ## tiles: one square per patient, 5 per row, coloured by TRG score
  tiles <- d[order(d$Group, d$TRG_score, d$SampleID), ]
  tiles$Index <- ave(seq_len(nrow(tiles)), tiles$Group, FUN = seq_along) - 1
  tiles$Group <- factor(paste0(carb, "-", tolower(tiles$Group), "\n", group_lab[as.character(tiles$Group)]),
                        levels = paste0(carb, "-", c("high", "low"), "\n", group_lab))
  waffle <- ggplot(tiles, aes(Index %% 5, -(Index %/% 5), fill = factor(TRG_score, levels = 0:3))) +
    geom_tile(width = 0.92, height = 0.92, colour = "grey40", linewidth = 0.3) +
    facet_wrap(~Group, ncol = 1) +
    coord_equal() +
    scale_fill_manual(values = trg_score_colors, name = "TRG score", drop = FALSE) +
    labs(caption = p_label(p, prefix = "p")) +
    theme_void(base_size = 9) +
    theme(strip.text = element_text(size = 7.5, lineheight = 1.1), plot.caption = element_text(hjust = 0.5, size = 9),
          legend.position = "right")
  save_panel(scatter + waffle + plot_layout(widths = c(2.4, 1)), name, width = 6.4, height = 3.4)
}
balance_panel("Sucrose",
              x_lab = "PWY-7238: sucrose biosynthesis II",
              y_lab = "PWY-621 + PWY-5384: sucrose degradation III + IV",
              group_lab = c(High = "PWY-7238 > 30000\n& PWY-621 + PWY-5384 = 0", Low = "PWY-7238 ≤ 30000\nOR PWY-621 + PWY-5384 > 0"),
              synthesis_cut = 30000, degradation_cut = 0, name = "SFig2B_sucrose_balance")
balance_panel("Rhamnose",
              x_lab = "DTDPRHAMSYN-PWY: dTDP-L-rhamnose biosynthesis",
              y_lab = "RHAMCAT-PWY: L-rhamnose degradation I",
              group_lab = c(High = "DTDPRHAMSYN-PWY ≥ 35000\nOR RHAMCAT-PWY ≤ 5000", Low = "DTDPRHAMSYN-PWY < 35000\n& RHAMCAT-PWY > 5000"),
              synthesis_cut = 35000, degradation_cut = 5000, name = "SFig2C_rhamnose_balance")


## ---- 4. Supplementary Figure 2D-M: taxon contributions ----------------------------

genus_colors <- c(
  Anaerostipes = "#A3C37B", Streptococcus = "#D96C6C", Escherichia = "#4C91E0", Bacteroides = "#F1A7C1",
  Klebsiella = "#D84B6A", Bifidobacterium = "#F0C16B", Roseburia = "#6F8A3F", Blautia = "#9A4E97",
  Faecalibacterium = "#8BC34A", Collinsella = "#1C85C4", Eubacterium = "#66B2E4", Dorea = "#FF9800",
  Enterococcus = "#66A3D2", Lactococcus = "#A0D468", Raoultella = "#E8C13B", Prevotella = "#78C38F",
  Lactobacillus = "#E17A6A", Others = "#BFBFBF", Unclassified = "#DDDDDD"
)
bacteroides_colors <- c(
  Bacteroides_thetaiotaomicron = "#FFD1DC", Bacteroides_xylanisolvens = "#D16A7D", Bacteroides_ovatus = "#FAD0C4",
  Bacteroides_uniformis = "#EC8B99", Bacteroides_caccae = "#FBCFE8", Bacteroides_faecis = "#A95C68",
  Bacteroides_fragilis = "#E56B6F", Bacteroides_cellulosilyticus = "#FFB3BA", Bacteroides_nordii = "#8E5B66",
  Bacteroides_stercoris = "#D97A8C", Bacteroides_intestinalis = "#E8A1B5", Bacteroides_salyersiae = "#C97B8E",
  Bacteroides_eggerthii = "#9F5968", Phocaeicola_dorei = "#C14450", Phocaeicola_plebeius = "#C9A9A9",
  Phocaeicola_vulgatus = "#F1A7C1", `Unresolved Bacteroides spp.` = "#DDDDDD", Others = "#BFBFBF"
)
blautia_colors <- c(
  Blautia_hydrogenotrophica = "#4F1D6B", Ruminococcus_gnavus = "#6A2E8C", Blautia_hansenii = "#8443A3",
  Blautia_wexlerae = "#9A4E97", Ruminococcus_torques = "#C16EB8", Blautia_obeum = "#E0A9DA",
  `Unresolved Blautia spp.` = "#DDDDDD", Others = "#BFBFBF"
)
## HUMAnN (older MetaPhlAn taxonomy) -> current MetaPhlAn species names
phocaeicola_alias <- c(Bacteroides_vulgatus = "Phocaeicola_vulgatus", Bacteroides_dorei = "Phocaeicola_dorei",
                       Bacteroides_plebeius = "Phocaeicola_plebeius")
harmonise <- function(x) ifelse(x %in% names(phocaeicola_alias), phocaeicola_alias[x], x)

target_ids <- c("PWY-621", "PWY-5384", "PWY-7238", "RHAMCAT-PWY", "DTDPRHAMSYN-PWY")
pathway_of_row <- sub(":.*$", "", sub("\\|.*$", "", rownames(abundance)))
strata <- which(!unstratified & pathway_of_row %in% target_ids)
contrib <- abundance[strata, baseline$SampleID]
contrib[coverage[strata, baseline$SampleID] < coverage_cut] <- 0
taxon <- sub("^[^|]*\\|", "", rownames(contrib))

long <- data.frame(
  Pathway = pathway_of_row[strata],
  Genus   = ifelse(grepl("g__", taxon), sub(".*g__([^.]*)\\..*", "\\1", taxon), "Unclassified"),
  Species = harmonise(ifelse(grepl("s__", taxon), sub(".*s__", "", taxon), "Unclassified")),
  contrib, check.names = FALSE
) %>%
  pivot_longer(all_of(baseline$SampleID), names_to = "SampleID", values_to = "Value") %>%
  left_join(baseline[, c("SampleID", "Response")], by = "SampleID")
group_n <- table(baseline$Response)

## Stacked-bar helper: group means (sum / group size) of the contributing taxa
contribution_bar <- function(grouped, palette, title, y_lab, legend_title, italic = TRUE) {
  grouped$Taxon <- factor(grouped$Taxon, levels = rev(intersect(names(palette), unique(grouped$Taxon))))
  ggplot(grouped, aes(Response, Value, fill = Taxon)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = palette, labels = function(x) gsub("_", " ", x), name = legend_title) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       labels = function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    labs(x = NULL, y = y_lab, title = title) +
    theme_tnt(base_size = 9, legend = "right") +
    theme(plot.title = element_text(size = 8.5),
          legend.text = element_text(size = 7, face = if (italic) "italic" else "plain"),
          legend.key.size = unit(0.35, "cm"))
}

## S2D-G: genus contributions (top-6 genera + Others + Unclassified)
genus_panels <- c("PWY-621" = "SFig2D", "PWY-5384" = "SFig2E", "RHAMCAT-PWY" = "SFig2F", "DTDPRHAMSYN-PWY" = "SFig2G")
pathway_titles <- setNames(sub("^[^:]*: ", "", rownames(raw)[match(names(genus_panels), pathway_ids)]), names(genus_panels))
for (id in names(genus_panels)) {
  part <- long %>% filter(Pathway == id, Value > 0)
  top_genera <- part %>% filter(Genus != "Unclassified") %>% group_by(Genus) %>%
    summarise(Value = sum(Value), .groups = "drop") %>% slice_max(Value, n = 6) %>% pull(Genus)   # <-- 6 genera shown
  grouped <- part %>%
    mutate(Taxon = case_when(Genus %in% top_genera ~ Genus, Genus == "Unclassified" ~ "Unclassified", TRUE ~ "Others")) %>%
    group_by(Response, Taxon) %>% summarise(Value = sum(Value), .groups = "drop") %>%
    mutate(Value = Value / as.numeric(group_n[as.character(Response)]))
  stopifnot(all(grouped$Taxon %in% names(genus_colors)))
  save_panel(contribution_bar(grouped, genus_colors, paste0(id, ": ", pathway_titles[[id]]),
                              "Functional abundance", "Genus"),
             paste0(genus_panels[[id]], "_", id, "_genus_contributions"), width = 3.4, height = 3.2)
}

## S2H-M: species contributions within Bacteroides / Blautia + species abundance
species <- read_metaphlan("species", samples = baseline$SampleID, drop_empty = FALSE)   # 26 x species (%)

for (genus in c("Bacteroides", "Blautia")) {
  if (genus == "Bacteroides") {
    in_genus <- grepl("^Bacteroides_|^Phocaeicola_(vulgatus|dorei|plebeius)$", colnames(species))
    palette <- bacteroides_colors; member_paths <- c("PWY-621", "RHAMCAT-PWY", "DTDPRHAMSYN-PWY")
    panel_names <- c("SFig2H", "SFig2I", "SFig2J")
  } else {
    in_genus <- grepl("^Blautia_|^Ruminococcus_(gnavus|torques)$", colnames(species))
    palette <- blautia_colors; member_paths <- c("PWY-7238", "RHAMCAT-PWY", "DTDPRHAMSYN-PWY")
    panel_names <- c("SFig2K", "SFig2L", "SFig2M")
  }
  part <- long %>% filter(Genus == genus, Value > 0) %>%
    mutate(Taxon = ifelse(Species %in% colnames(species)[in_genus], Species, paste0("Unresolved ", genus, " spp.")))
  functional_species <- unique(part$Taxon[part$Pathway %in% member_paths])

  for (j in 1:2) {
    id <- c("RHAMCAT-PWY", "DTDPRHAMSYN-PWY")[j]
    grouped <- part %>% filter(Pathway == id) %>% group_by(Response, Taxon) %>%
      summarise(Value = sum(Value), .groups = "drop") %>%
      mutate(Value = Value / as.numeric(group_n[as.character(Response)]))
    stopifnot(all(grouped$Taxon %in% names(palette)))
    save_panel(contribution_bar(grouped, palette, id, "Functional abundance", "Species"),
               paste0(panel_names[j], "_", genus, "_", id), width = 3.6, height = 3.2)
  }
  ## species relative abundance of the same taxa (MetaPhlAn)
  sp_long <- data.frame(SampleID = rownames(species), species[, in_genus], check.names = FALSE) %>%
    pivot_longer(-SampleID, names_to = "Species", values_to = "Value") %>%
    left_join(baseline[, c("SampleID", "Response")], by = "SampleID") %>%
    mutate(Taxon = ifelse(Species %in% functional_species, Species, "Others")) %>%
    group_by(Response, Taxon) %>% summarise(Value = sum(Value), .groups = "drop") %>%
    mutate(Value = Value / as.numeric(group_n[as.character(Response)]))
  stopifnot(all(sp_long$Taxon %in% names(palette)))
  save_panel(contribution_bar(sp_long, palette, "Species abundance", "Mean relative abundance (%)", "Species"),
             paste0(panel_names[3], "_", genus, "_species_abundance"), width = 3.6, height = 3.2)
}

message("Done: 06b_SFig2_MetaCyc_pathways.R")
