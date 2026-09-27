## =============================================================================
##  09c_SFig4_Random_forest_ROC_SHAP.R
##  ---------------------------------------------------------------------------
##  PANELS      Supplementary Figure 4B    ROC curves of the nine random-forest models
##              Supplementary Figure 4C-J  top-10 features per model: mean |SHAP|
##                                         bars + MDA / MDG / consistency / score tiles
##
##  PURPOSE     Display the results computed (or replayed) by script 09b.  No
##              model is fitted here.
##
##  INPUT       output/cache/random_forest_results.rds   (from 09b)
##
##  OUTPUT      output/figures/SFig4B_random_forest_ROC.{svg,png}
##              output/figures/SFig4C_SGB_SHAP ... SFig4J_KO_SHAP  (.svg/.png)
##              output/figures/SFig4_SHAP_tile_legend.{svg,png}
##              output/tables/SFig4CJ_top10_features.csv
##
##  METHODS     ROC curves from pooled out-of-fold P(pCR) (pROC 1.19).  Feature
##              panels show the ten features with the highest mean |SHAP|; the
##              tiles are the fold-averaged percentile ranks (0-1) of mean
##              decrease in accuracy (MDA) and Gini (MDG) - NOT percentages,
##              although the submitted axis label says "(%)" - the screening
##              consistency (0.7-1) and -log10 of the screening score (capped
##              at 15).  The score is a screening quantity, not a calibrated P.
##
##  R PACKAGES  pROC, ggplot2, patchwork
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(pROC)
  library(patchwork)
})
report_versions(c("pROC", "ggplot2", "patchwork"))

all_results <- readRDS(cache_file("random_forest_results.rds"))

## Display names and colours of the nine models (as submitted)
input_labels <- c(Strain = "SGB", Species = "Species", Genus = "Genus", Family = "Family", Phylum = "Phylum",
                  Metabolite = "Metabolite", Pathway = "MetaCyc pathway", KO = "KEGG ortholog", CI = "cStage + CEA")
input_colors <- c(Strain = "#3E6DB5", Species = "#2CB1C9", Genus = "#5FC8B6", Family = "#98DDD9", Phylum = "#C7E6F2",
                  Metabolite = "#5FBF66", Pathway = "#F2B233", KO = "#D48806", CI = "#AEAEAE")
model_order <- names(input_labels)


## ---- 1. Supplementary Figure 4B: ROC curves ---------------------------------------------

roc_df <- do.call(rbind, lapply(model_order, function(nm) {
  r <- all_results[[nm]]$ROC
  data.frame(Model = nm, FPR = 1 - r$specificities, TPR = r$sensitivities)[order(1 - r$specificities, r$sensitivities), ]
}))
legend_labels <- sapply(model_order, function(nm) sprintf("%s (AUC = %.3f)", input_labels[[nm]], all_results[[nm]]$AUC))
roc_df$Model <- factor(roc_df$Model, levels = model_order)

p_roc <- ggplot(roc_df, aes(FPR, TPR, colour = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.4) +
  geom_path(linewidth = 0.9, alpha = 0.9) +
  scale_colour_manual(values = input_colors, labels = legend_labels, name = "Input data") +
  scale_x_continuous(expand = c(0.01, 0)) + scale_y_continuous(expand = c(0.01, 0)) +
  coord_equal() +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_tnt(base_size = 10, legend = "right") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        legend.text = element_text(size = 8), legend.key.width = unit(0.8, "cm"))
save_panel(p_roc, "SFig4B_random_forest_ROC", width = 6.4, height = 4)


## ---- 2. Supplementary Figure 4C-J: SHAP bars + importance tiles --------------------------

panel_letters <- c(Strain = "C", Species = "D", Genus = "E", Family = "F", Phylum = "G",
                   Metabolite = "H", Pathway = "I", KO = "J")
panel_titles  <- c(input_labels[names(panel_letters)]); panel_titles[c("Pathway", "KO")] <- c("MetaCyc pathways", "KEGG orthologs")

## Submitted display names for KOs whose HUMAnN name is missing / generic
ko_labels <- c(K19701 = "aminopeptidase YwaD", K19140 = "CRISPR-associated protein Csm5",
               K06926 = "uncharacterized protein", K01932 = "gamma-polyglutamate synthase",
               K06919 = "putative DNA primase/helicase")

## Tile columns: variable, axis label, colour, display limits and legend breaks
tiles_spec <- list(
  MDA   = list(label = "MDA (%)",     colour = "#2C7FB8", limits = c(0, 1),    breaks = c(0, 0.25, 0.5, 0.75, 1)),
  MDG   = list(label = "MDG (%)",     colour = "#2F9EAA", limits = c(0, 1),    breaks = c(0, 0.25, 0.5, 0.75, 1)),
  Consistency = list(label = "Consistency", colour = "#4C956C", limits = c(0.7, 1), breaks = c(0.7, 0.8, 0.9, 1)),
  Score = list(label = "-Log10(p)",   colour = "#D6406B", limits = c(0, 15),   breaks = c(0, 5, 10, 15))
)

feature_label <- function(nm, features) {
  if (nm %in% c("Strain", "Species", "Genus")) return(taxon_plotmath(features))
  lab <- gsub("_", " ", features)
  if (nm == "KO") {
    ids <- sub(":.*$", "", features)
    lab <- sub(" \\[EC:.*$", "", lab)
    hit <- ids %in% names(ko_labels)
    lab[hit] <- paste0(ids[hit], ": ", ko_labels[ids[hit]])
  }
  sprintf('"%s"', lab)                                          # plain plotmath strings
}

shap_panel <- function(nm) {
  d <- all_results[[nm]]$Feature_summary
  d <- head(d[order(-d$MeanAbsSHAP), ], 10)                    # <-- number of features shown
  d$Score <- -log10(pmax(d$Pvalue, .Machine$double.xmin))
  d$Feature <- factor(d$Feature, levels = rev(d$Feature))
  labs_pm <- setNames(feature_label(nm, levels(d$Feature)), levels(d$Feature))

  bar <- ggplot(d, aes(MeanAbsSHAP, Feature)) +
    geom_col(fill = input_colors[[nm]], width = 0.75) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    scale_y_discrete(labels = function(b) parse(text = labs_pm[b])) +
    labs(x = "Mean |SHAP|", y = NULL, title = panel_titles[[nm]]) +
    theme_tnt(base_size = 9) +
    theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),      # light vertical grid as submitted
          axis.ticks.y = element_blank(), plot.title = element_text(size = 10))

  tiles <- lapply(names(tiles_spec), function(v) {
    sp <- tiles_spec[[v]]
    ggplot(d, aes(1, Feature, fill = .data[[v]])) +
      geom_tile(colour = NA) +
      scale_fill_gradient(low = "white", high = sp$colour, limits = sp$limits, oob = scales::squish, guide = "none") +
      scale_x_continuous(expand = c(0, 0)) +
      labs(x = sp$label, y = NULL) +
      theme_void(base_size = 8) +
      theme(axis.title.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7.5),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6))
  })
  wrap_plots(c(list(bar), tiles), nrow = 1, widths = c(6, rep(0.55, 4)))
}

top10 <- list()
for (nm in names(panel_letters)) {
  n_feat <- min(10, nrow(all_results[[nm]]$Feature_summary))
  save_panel(shap_panel(nm), paste0("SFig4", panel_letters[[nm]], "_", nm, "_SHAP"),
             width = 5.2, height = 1.3 + 0.27 * n_feat)
  top10[[nm]] <- cbind(Input = nm, head(all_results[[nm]]$Feature_summary[order(-all_results[[nm]]$Feature_summary$MeanAbsSHAP), ], 10))
}
save_table(do.call(rbind, top10), "SFig4CJ_top10_features")

## Shared colour-bar legend for the four tile columns (one small panel)
legend_plots <- lapply(names(tiles_spec), function(v) {
  sp <- tiles_spec[[v]]
  ggplot(data.frame(x = seq(sp$limits[1], sp$limits[2], length.out = 100), y = 1), aes(x, y, fill = x)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = sp$colour, limits = sp$limits, guide = "none") +
    scale_x_continuous(breaks = sp$breaks, expand = c(0, 0)) +
    labs(x = NULL, y = sp$label) +
    theme_tnt(base_size = 8) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4))
})
save_panel(wrap_plots(legend_plots, ncol = 1), "SFig4_SHAP_tile_legend", width = 2.8, height = 2.4)

message("Done: 09c_SFig4_Random_forest_ROC_SHAP.R")
