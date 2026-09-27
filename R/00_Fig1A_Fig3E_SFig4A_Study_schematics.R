## =============================================================================
##  00_Fig1A_Fig3E_SFig4A_Study_schematics.R
##  ---------------------------------------------------------------------------
##  PANELS      Figure 1A (data-availability intersection plot + study timeline),
##              Figure 3E (tryptophan catabolism scheme),
##              Supplementary Figure 4A (random-forest workflow scheme)
##
##  PURPOSE     Reproduce the three "schematic" panels.  The intersection
##              (UpSet-style) plot of Figure 1A is fully data-driven: it counts,
##              for each of the 26 patients, which omic layers are available at
##              baseline and after radiotherapy (RT).  The treatment timeline,
##              the tryptophan pathway and the model workflow are diagrams that
##              were finished in Adobe Illustrator for the submission; the R
##              versions here contain the same text, counts and topology but
##              use simple boxes instead of the hand-drawn icons and chemical
##              structures.
##
##  INPUT       input/metadata/stool_samples.csv     42 stool metagenomes
##              input/metadata/rna_samples.csv       27 tumor transcriptomes
##              input/metabolome/fecal_metabolites.csv  32 fecal metabolomes
##              input/metadata/patients.csv          response labels
##
##  OUTPUT      output/figures/Fig1A_data_availability_upset.{svg,png}
##              output/figures/Fig1A_study_timeline.{svg,png}
##              output/figures/Fig3E_tryptophan_pathway.{svg,png}
##              output/figures/SFig4A_random_forest_workflow.{svg,png}
##              output/tables/Fig1A_data_availability.csv
##
##  METHODS     No statistics.  Set intersections are computed with base R.
##
##  R PACKAGES  ggplot2, patchwork, svglite
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages(library(patchwork))
report_versions(c("ggplot2", "patchwork", "svglite"))


## ---- 1. Figure 1A: which omic layers exist for each patient? ----------------

patients <- read_patients()
stool    <- read_stool_samples()
metab    <- read_metabolites()
rna      <- read.csv(input_file("metadata", "rna_samples.csv"), stringsAsFactors = FALSE)

## One logical column per (layer, time point); rows are patients.
layer_sets <- list(
  "Metagenome|Baseline"  = stool$SubjectID[stool$Timepoint == "Before"],
  "Metagenome|After RT"  = stool$SubjectID[stool$Timepoint == "Ongoing"],
  "Metabolome|Baseline"  = metab$SubjectID[metab$Timepoint == "Before"],
  "Metabolome|After RT"  = metab$SubjectID[metab$Timepoint == "Ongoing"],
  "Host RNA-seq|Baseline" = rna$SubjectID[rna$Timepoint == "Before"],
  "Host RNA-seq|After RT" = rna$SubjectID[rna$Timepoint == "Ongoing"]
)
membership <- sapply(layer_sets, function(ids) patients$SubjectID %in% ids)
rownames(membership) <- patients$SubjectID
set_sizes <- colSums(membership)

cat("Omic profiles:", sum(set_sizes), "=",
    sum(set_sizes[1:2]), "metagenomes +", sum(set_sizes[3:4]), "metabolomes +",
    sum(set_sizes[5:6]), "transcriptomes\n")
stopifnot(sum(set_sizes) == 101)

## Intersections: identical availability patterns are grouped and counted.
pattern <- apply(membership, 1, function(z) paste(as.integer(z), collapse = ""))
combos  <- as.data.frame(table(pattern), stringsAsFactors = FALSE)
names(combos) <- c("pattern", "n_patients")
combos$degree <- nchar(gsub("0", "", combos$pattern))

## Display order: many-layer patterns first, then by frequency  <-- change here
## to reorder the columns of the intersection plot.
combos <- combos[order(-combos$degree, -combos$n_patients, combos$pattern), ]
combos$x <- seq_len(nrow(combos))

matrix_long <- do.call(rbind, lapply(seq_len(nrow(combos)), function(i) {
  data.frame(x = combos$x[i], set = names(layer_sets),
             member = strsplit(combos$pattern[i], "")[[1]] == "1")
}))
matrix_long$set <- factor(matrix_long$set, levels = rev(names(layer_sets)))
segments <- do.call(rbind, lapply(split(matrix_long[matrix_long$member, ], matrix_long$x[matrix_long$member]),
  function(d) data.frame(x = d$x[1], ymin = min(as.integer(d$set)), ymax = max(as.integer(d$set)))))

save_table(cbind(SubjectID = rownames(membership), as.data.frame(membership)), "Fig1A_data_availability")

## Design: bars and active dots in slate grey (#6E7C89), inactive dots in light grey.
dot_on  <- "#6E7C89"
dot_off <- "#D7DDE1"

p_top <- ggplot(combos, aes(x, n_patients)) +
  geom_col(fill = dot_on, width = 0.6) +
  geom_text(aes(label = n_patients), vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_x_continuous(expand = expansion(add = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25)), breaks = c(0, 2, 4, 6)) +
  labs(x = NULL, y = "Available samples (n)") +
  theme_tnt(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())

p_matrix <- ggplot(matrix_long, aes(x, set)) +
  geom_segment(data = segments, aes(x = x, xend = x, y = ymin, yend = ymax), inherit.aes = FALSE,
               colour = dot_on, linewidth = 0.9) +
  geom_point(aes(colour = member), size = 3.2) +
  scale_colour_manual(values = c("TRUE" = dot_on, "FALSE" = dot_off), guide = "none") +
  scale_x_continuous(expand = expansion(add = 0.6)) +
  scale_y_discrete(labels = function(z) sub(".*\\|", "", z)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 9) +
  theme(axis.text.y = element_text(size = 7.5, hjust = 1))

## layer names (Metagenome / Metabolome / Host RNA-seq) bracketing two rows each,
## drawn as a narrow plot that shares the y scale of the dot matrix
layer_labels <- data.frame(
  set   = factor(names(layer_sets), levels = rev(names(layer_sets))),
  layer = sub("\\|.*", "", names(layer_sets))
)
bracket <- data.frame(layer = unique(layer_labels$layer), y = c(5.5, 3.5, 1.5))   # rows are 6..1 from the top
p_left <- ggplot(layer_labels, aes(x = 1, y = set)) +
  geom_blank() +
  geom_segment(data = bracket, aes(x = 1, xend = 1, y = y - 0.45, yend = y + 0.45), inherit.aes = FALSE,
               colour = "grey40", linewidth = 0.4) +
  geom_text(data = bracket, aes(x = 0.95, y = y, label = layer), inherit.aes = FALSE, hjust = 1, size = 3) +
  scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  theme_void()

p_side <- ggplot(data.frame(set = factor(names(set_sizes), levels = rev(names(layer_sets))), n = set_sizes),
                 aes(n, set)) +
  geom_col(fill = "#AEBED0", width = 0.6) +
  geom_text(aes(label = n), hjust = -0.3, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "Patients with data (n)", y = NULL) +
  theme_tnt(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

p_upset <- (plot_spacer() + p_top + plot_spacer()) / (p_left + p_matrix + p_side) +
  plot_layout(widths = c(0.9, 3, 1.1), heights = c(1.1, 1))
save_panel(p_upset, "Fig1A_data_availability_upset", width = 6.8, height = 3.6)


## ---- 2. Figure 1A: treatment timeline (text version of the artwork) ----------

n_resp <- table(patients$Response)
timeline <- data.frame(
  x = c(0, 1, 2, 3), y = 1,
  label = c("LARC patients\n(N = 26)", "RT\n(25 Gy / 5 fx)", "CAPOX\n(4 cycles)", "Surgery\n(TME)")
)
p_timeline <- ggplot(timeline, aes(x, y)) +
  annotate("segment", x = -0.1, xend = 3.1, y = 1, yend = 1, linewidth = 1.2, colour = "grey30") +
  geom_label(aes(label = label), size = 3.2, label.padding = unit(0.35, "lines"), fill = "white") +
  annotate("text", x = 0.5, y = 0.7, label = "Baseline\nsampling", size = 2.8, colour = "grey30") +
  annotate("text", x = 1.5, y = 0.7, label = "After-RT\nsampling", size = 2.8, colour = "grey30") +
  annotate("text", x = 3.6, y = 1.12, hjust = 0, size = 3.2, colour = response_colors[["pCR"]],
           label = sprintf("TRG = 0  ->  pCR (n = %d)", n_resp[["pCR"]])) +
  annotate("text", x = 3.6, y = 0.88, hjust = 0, size = 3.2, colour = response_colors[["non-pCR"]],
           label = sprintf("TRG > 0  ->  non-pCR (n = %d)", n_resp[["non-pCR"]])) +
  annotate("text", x = 1.5, y = 1.45, label = "Total neoadjuvant therapy (TNT)", fontface = "bold", size = 3.8) +
  coord_cartesian(xlim = c(-0.4, 5.6), ylim = c(0.5, 1.6), clip = "off") +
  theme_void()
save_panel(p_timeline, "Fig1A_study_timeline", width = 7.5, height = 2.2)


## ---- 3. Figure 3E: microbial tryptophan catabolism through indole-3-pyruvate --
## Reductive branch (left, teal) yields ILA and IPA; decarboxylative branch
## (right, salmon) yields IAA and skatole.  Chemical structures are not drawn.

nodes <- data.frame(
  x = c(1, 1, 0, 2, 0, 2, 0, 2),
  y = c(6, 5, 4, 4, 3, 3, 2, 2),
  label = c("Tryptophan", "Indole-3-pyruvic acid",
            "Indole-3-lactic acid (ILA)", "Indole-3-acetaldehyde",
            "Indole-3-acrylic acid", "Indole-3-acetic acid (IAA)",
            "Indole-3-propionic acid (IPA)", "Skatole"),
  bold = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE)
)
arrows <- data.frame(from = c(1, 2, 2, 3, 4, 5, 6), to = c(2, 3, 4, 5, 6, 7, 8))
arrows$x    <- nodes$x[arrows$from];  arrows$xend <- nodes$x[arrows$to]
arrows$y    <- nodes$y[arrows$from] - 0.22; arrows$yend <- nodes$y[arrows$to] + 0.22
p_trp <- ggplot() +
  geom_segment(data = arrows, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.5) +
  geom_text(data = nodes, aes(x, y, label = label, fontface = ifelse(bold, "bold", "plain")), size = 3.4) +
  annotate("text", x = 0, y = 4.6, label = "Reductive", colour = response_colors[["pCR"]], fontface = "bold", size = 4) +
  annotate("text", x = 2, y = 4.6, label = "Decarboxylative", colour = response_colors[["non-pCR"]], fontface = "bold", size = 4) +
  coord_cartesian(xlim = c(-0.8, 2.8), ylim = c(1.6, 6.3), clip = "off") +
  theme_void()
save_panel(p_trp, "Fig3E_tryptophan_pathway", width = 5.5, height = 4.2)


## ---- 4. Supplementary Figure 4A: random-forest workflow -------------------------
## Text boxes summarising the screening and cross-validation procedure that is
## implemented in scripts 09a-09c (see those headers for parameters).

steps <- data.frame(y = 5:1, label = c(
  "Baseline and after-RT profiles of 26 patients\n(pCR n = 11, non-pCR n = 15)",
  "Feature screening on the full cohort:\nconcordant pCR vs non-pCR direction at both time points",
  "Grid search over screening cut-offs\n(prevalence, effect size D, consistency, screening score, Wilcoxon P, top-N)",
  "Random forest (1,000 trees), stratified 5-fold CV grouped by patient;\nout-of-fold probabilities pooled for the ROC / AUC",
  "Feature importance: mean |SHAP| (100 simulations),\nfold-wise MDA / MDG percentile ranks"
))
p_rf <- ggplot(steps, aes(0, y)) +
  geom_segment(data = steps[1:4, ], aes(x = 0, xend = 0, y = y - 0.32, yend = y - 0.68),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.5) +
  geom_label(aes(label = label), size = 3.1, fill = "#F1F5F6", label.padding = unit(0.4, "lines"), linewidth = 0.3) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(0.5, 5.5), clip = "off") +
  theme_void()
save_panel(p_rf, "SFig4A_random_forest_workflow", width = 7, height = 5.5)

message("Done: 00_Fig1A_Fig3E_SFig4A_Study_schematics.R")
