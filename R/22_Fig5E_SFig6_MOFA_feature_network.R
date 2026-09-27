## =============================================================================
##  22_Fig5E_SFig6_MOFA_feature_network.R
##  ---------------------------------------------------------------------------
##  PANELS      Figure 5E              leading features of MOFA Factor 7 (lollipop
##                                     of normalised loadings, grouped by data type)
##              Supplementary Fig. 6   network of the same features connected by
##                                     age/BMI/sex-adjusted partial Spearman
##                                     correlations (nominal permutation P < 0.10)
##
##  PURPOSE     Show which species, KEGG orthologs and metabolites carry the
##              response-associated Factor 7, and whether those features
##              co-vary across patients independently of the factor model.
##
##  INPUT       output/cache/mofa_model.RData          (script 20: mofa_input, mofa_model)
##              output/cache/mofa_associations.RData   (script 21: a = factor associations)
##              output/cache/four_omics_inputs.RData   (script 18: KO names for labels)
##              input/metadata/patients_covariates.csv   Sex, Age, BMI (public; the
##                  restricted patients_clinical.csv is used instead when present)
##                  - without either file Figure 5E is still drawn; the adjusted
##                    network (Supplementary Figure 6) is skipped with a message.
##
##  OUTPUT      output/figures/Fig5E_Factor7_loadings, SFig6_Factor7_network  (.svg/.png)
##              output/tables/Fig5E_Factor7_features.csv
##              output/tables/SFig6_Factor7_partial_correlations.csv
##
##  METHODS     Feature selection (candidate_features, R/_helpers_mofa.R): for
##              Factor 7 a view must explain >= 1 % of its variance; a feature
##              must have |weight| >= 25 % of the largest weight of its view and
##              a factor-wide normalised |weight| > 0.10; at most five features
##              per view.  Loadings are divided by the largest |weight| of the
##              factor and oriented so that the pCR direction is positive.
##              Hedges' g of each feature (patient means, pCR - non-pCR) is
##              tabulated.  Network edges: partial Spearman correlation of the
##              patient means of two features (visits where both were observed)
##              after removing Age, BMI and Sex by rank regression; P from 999
##              permutations of the residuals (seed 20260908 + 300 + pair
##              index); BH q over the full candidate family.  The drawn graph is
##              the induced subgraph of features with normalised |loading|
##              > 0.20 (12 features in the frozen model; glutamic acid, with a
##              loading of 0.11-0.20, is tested but not drawn).
##
##  WHAT CHANGES WHAT
##              display_loading_cut (0.20)   which features are drawn - lowering
##                                           it adds glutamic acid to both panels
##              edge_P_cut (0.10)            which correlations become edges in SFig 6
##              mofa_cfg$minimum_* (helpers) the candidate family itself (and
##                                           therefore the BH q values)
##              layout seed                  only the arrangement of SFig 6 nodes
##
##  R PACKAGES  igraph, ggraph, ggplot2
## =============================================================================

source("R/_common.R")
source("R/_helpers_mofa.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(ggraph)
})
report_versions(c("dplyr", "igraph", "ggraph", "ggplot2"))

load(cache_file("mofa_model.RData"))                                  # mofa_input, mofa_model
load(cache_file("mofa_associations.RData"))                           # a, effects, ve
load(cache_file("four_omics_inputs.RData"))                           # coherence (KO names)

lead_factor        <- a$axes[1]                                       # "Factor7" in the frozen model
display_loading_cut <- 0.20                                           # <-- nodes drawn (normalised |loading|)
edge_P_cut         <- 0.10                                            # <-- edges drawn (nominal permutation P)

view_labels <- c(species = "Species", ko = "KEGG ortholog", metabolite = "Metabolite", host = "Host gene")
view_colors <- c(species = "#74ADA0", ko = "#CD9974", metabolite = "#7AA0BD", host = "#A492BD")
view_shapes <- c(species = 21, ko = 24, metabolite = 22, host = 23)


## ---- 1. Candidate features of the leading factor -----------------------------------------

nodes <- candidate_features(mofa_model$extracted$weights, mofa_model$extracted$variance, lead_factor)

## KO names for a refit without an annotation column (frozen weights already carry them)
ko_i <- which(nodes$view == "ko" & nodes$feature_label == nodes$feature)
if (length(ko_i)) nodes$feature_label[ko_i] <- display_feature_label("ko", nodes$feature[ko_i],
                                                                      coherence$ko_annotation$KO_name[match(nodes$feature[ko_i], coherence$ko_annotation$KO)])

## observed effect of every candidate: Hedges' g of the patient means (pCR - non-pCR)
md <- mofa_input$metadata
nodes$g <- vapply(seq_len(nrow(nodes)), function(i) {
  x <- mofa_input$values[[nodes$view[i]]][, nodes$feature[i]]
  u <- patient_means(x, md)
  hedges(u[, 1], md$Response[match(rownames(u), md$SubjectID)])
}, numeric(1))
nodes$oriented_loading <- nodes$factor_scaled_weight * factor_sign(a, lead_factor)   # pCR direction positive

cat(sprintf("%s: %d candidate features (%s)\n", lead_factor, nrow(nodes),
            paste(sprintf("%s %d", names(table(nodes$view)), table(nodes$view)), collapse = ", ")))


## ---- 2. Adjusted partial correlations within the candidate family ---------------------------

## Every pair of candidates is tested (the BH family), even those not drawn later.
## Each pair uses the visits at which both features were observed, averaged per patient.
have_clinical <- has_covariates(c("Age", "BMI", "Sex"))      # patients_covariates.csv (public) or patients_clinical.csv
edges <- NULL
if (have_clinical) {
  covariates <- read_patients(clinical = TRUE)[, c("SubjectID", "Age", "BMI", "Sex")]
  pairs <- combn(seq_len(nrow(nodes)), 2)
  edges <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
    i <- pairs[1, j]; k <- pairs[2, j]
    x <- mofa_input$values[[nodes$view[i]]][, nodes$feature[i]]
    y <- mofa_input$values[[nodes$view[k]]][, nodes$feature[k]]
    keep <- is.finite(x) & is.finite(y)
    u <- patient_means(cbind(x = x[keep], y = y[keep]), md[keep, ])
    cv <- covariates[match(rownames(u), covariates$SubjectID), ]
    r <- partial_spearman(u[, 1], u[, 2], cv, seed = mofa_cfg$seed + 300L + j)
    data.frame(from = nodes$key[i], to = nodes$key[k], from_label = nodes$feature_label[i], to_label = nodes$feature_label[k],
               n = r$n, rho = r$rho, P = r$P, t_P = r$t_P, df = r$df)
  }))
  edges$q <- p.adjust(edges$P, "BH")
  save_table(edges, "SFig6_Factor7_partial_correlations")
} else {
  message("Age / BMI / Sex not available (input/metadata/patients_covariates.csv): the adjusted network (SFig 6) is skipped.")
}


## ---- 3. Displayed features ---------------------------------------------------------------

shown <- nodes[abs(nodes$factor_scaled_weight) > display_loading_cut, ]
shown$view <- factor(shown$view, levels = view_order)
shown <- shown[order(shown$view, -shown$oriented_loading), ]
print(shown[, c("view", "feature_label", "weight", "oriented_loading", "g")], row.names = FALSE, digits = 3)
save_table(nodes[, c("view", "feature", "feature_label", "weight", "relative_weight", "factor_scaled_weight", "oriented_loading", "g", "r2")],
           "Fig5E_Factor7_features")

## plotmath labels: Latin binomials in italics, everything else upright
shown$label_pm <- ifelse(shown$view == "species", taxon_plotmath(shown$feature), sprintf('"%s"', shown$feature_label))


## ---- 4. Figure 5E: lollipop of normalised loadings, grouped by data type -----------------------

## y positions: one row per feature, a one-row gap between data types (top = species)
gap <- 1
shown$y <- NA_real_
pos <- 0
for (v in rev(levels(droplevels(shown$view)))) {                     # bottom-most view first
  i <- rev(which(shown$view == v))                                   # lowest loading at the bottom
  shown$y[i] <- pos + seq_along(i); pos <- pos + length(i) + gap
}
strips <- shown %>% dplyr::group_by(view) %>% dplyr::summarise(ymin = min(y) - 0.4, ymax = max(y) + 0.4, n = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(label = view_labels[as.character(view)], text_size = ifelse(n >= 4, 3, 2.2))   # short strips get a smaller label
x_max <- max(1, max(abs(shown$oriented_loading)))
strip_x <- c(x_max * 1.08, x_max * 1.16)

p_5e <- ggplot(shown, aes(oriented_loading, y)) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 1.2) +
  geom_segment(aes(x = 0, xend = oriented_loading, yend = y, colour = view), linewidth = 0.7) +
  geom_point(aes(fill = view, shape = view), size = 3.6, colour = "grey25", stroke = 0.5) +
  annotate("rect", xmin = strip_x[1], xmax = strip_x[2], ymin = strips$ymin, ymax = strips$ymax, fill = view_colors[as.character(strips$view)]) +
  annotate("text", x = mean(strip_x), y = (strips$ymin + strips$ymax) / 2, label = strips$label, angle = 270, colour = "white", size = strips$text_size) +
  scale_colour_manual(values = view_colors, guide = "none") +
  scale_fill_manual(values = view_colors, guide = "none") +
  scale_shape_manual(values = view_shapes, guide = "none") +
  scale_y_continuous(breaks = shown$y, labels = function(b) parse(text = shown$label_pm[match(b, shown$y)]), expand = expansion(add = 0.6)) +
  scale_x_continuous(breaks = seq(-1, 1, 0.5), expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(xlim = c(min(0, min(shown$oriented_loading)) - 0.02, strip_x[2]), clip = "off") +
  labs(x = "Normalized feature loading", y = NULL, title = sub("Factor", "Factor ", lead_factor)) +
  theme_tnt(base_size = 10) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 11))
save_panel(p_5e, "Fig5E_Factor7_loadings", width = 4.8, height = 4)


## ---- 5. Supplementary Figure 6: adjusted feature network ----------------------------------------

if (have_clinical) {
  draw <- edges[edges$from %in% shown$key & edges$to %in% shown$key & is.finite(edges$P) & edges$P < edge_P_cut, ]
  cat(sprintf("SFig 6: %d of %d within-display pairs with P < %.2f\n", nrow(draw),
              sum(edges$from %in% shown$key & edges$to %in% shown$key), edge_P_cut))
  print(draw[order(draw$P), c("from_label", "to_label", "n", "rho", "P", "q")], row.names = FALSE, digits = 3)

  vertices <- data.frame(name = shown$key, view = as.character(shown$view), label_pm = shown$label_pm,
                         loading = abs(shown$oriented_loading), stringsAsFactors = FALSE)
  graph <- graph_from_data_frame(draw[, c("from", "to", "rho", "P")], directed = FALSE, vertices = vertices)

  set.seed(mofa_cfg$seed)                                            # <-- layout only (the submitted panel was arranged by hand)
  lay <- create_layout(graph, layout = "fr", weights = abs(E(graph)$rho))

  p_s6 <- ggraph(lay) +
    geom_edge_link(aes(edge_width = abs(rho)), colour = "grey55", alpha = 0.8) +
    scale_edge_width_continuous(range = c(0.5, 3), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75), name = "Partial Spearman") +
    geom_node_point(aes(fill = view, shape = view), size = 6, colour = "grey20", stroke = 0.6) +
    geom_node_text(aes(label = label_pm), parse = TRUE, repel = TRUE, seed = mofa_cfg$seed, size = 3.1,
                   max.overlaps = Inf, box.padding = 0.9, point.padding = 0.6, force = 3, max.time = 5, max.iter = 10000,
                   segment.colour = "grey30", segment.size = 0.3) +
    scale_fill_manual(values = view_colors, labels = view_labels, breaks = c("species", "metabolite", "ko", "host"), name = "Data type") +
    scale_shape_manual(values = view_shapes, labels = view_labels, breaks = c("species", "metabolite", "ko", "host"), name = "Data type") +
    guides(fill = guide_legend(override.aes = list(size = 3), order = 1), shape = guide_legend(order = 1),
           edge_width = guide_legend(order = 2)) +
    coord_cartesian(clip = "off") +
    theme_void(base_size = 10) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.direction = "horizontal",
          legend.title.position = "left", legend.title = element_text(size = 10, margin = margin(r = 8)),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = 0.3),
          legend.margin = margin(3, 8, 3, 8), legend.box.just = "left", plot.margin = margin(15, 25, 5, 25))
  save_panel(p_s6, "SFig6_Factor7_network", width = 7.5, height = 6.5)
}

message("Done: 22_Fig5E_SFig6_MOFA_feature_network.R")
