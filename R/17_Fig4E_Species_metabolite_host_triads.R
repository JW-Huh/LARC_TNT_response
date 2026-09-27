## =============================================================================
##  17_Fig4E_Species_metabolite_host_triads.R
##  ---------------------------------------------------------------------------
##  PANEL       Figure 4E  six species - metabolite - host gene "triads": the three
##                         pairwise Spearman correlations drawn as a triangle
##
##  PURPOSE     Illustrate concrete three-layer associations (a species, a
##              metabolite it may produce/consume, and a tumour gene) that were
##              picked from the correlation analyses of Figures 3C, 4C and 4D.
##              The six triads are a fixed display list from the manuscript; no
##              ranking statistic or FDR is recomputed here.
##
##  INPUT       output/cache/host_metabolite_inputs.RData  (script 14: 17 patients)
##              input/microbiome/metaphlan_species.csv
##              input/metadata/legacy_paired_subjects.csv   the 10 patients whose
##                  species profiles entered the original triad analysis
##
##  OUTPUT      output/figures/Fig4E_triads.{svg,png}
##              output/tables/Fig4E_triad_edges.csv, Fig4E_triad_nodes.csv
##
##  METHODS     Spearman rho / P (normal approximation) for each pair.  NOTE the
##              sample sizes differ by edge type, as in the submitted figure:
##              metabolite-gene edges use the 17 matched patients, whereas the
##              two species edges use only the 10 patients of the legacy paired
##              subset that also have a biopsy (n = 7-10).  Node colour = log2
##              ratio pCR / non-pCR (species: means; metabolite: medians; gene:
##              DESeq2 log2FC).  Solid line = P < 0.05, dashed = not.  These are
##              associations, not a causal chain or mediation analysis.
##
##  R PACKAGES  ggplot2
## =============================================================================

source("R/_common.R")
report_versions("ggplot2")

load(cache_file("host_metabolite_inputs.RData"))                      # meta, met_raw, met_log, host_expr, host_deg
legacy <- read_csv_na(input_file("metadata", "legacy_paired_subjects.csv"))$SubjectID

## species (raw %) for the legacy subset of patients that also have a biopsy
sp_meta <- meta[meta$SubjectID %in% legacy, ]
species_raw <- read_metaphlan("species", samples = sp_meta$SampleID, drop_empty = FALSE)
rownames(species_raw) <- sp_meta$RNA_sample_id
species_log <- log10(species_raw + 1e-6)
cat("Species edges use", nrow(species_raw), "patients; metabolite-gene edges use", nrow(meta), "\n")

## The six displayed triads  <-- edit to draw other combinations
triads <- data.frame(
  Species    = c("Coprococcus_comes", "Coprococcus_comes", "Coprococcus_comes", "Roseburia_inulinivorans", "Ruminococcus_lactaris", "Bacteroides_fragilis"),
  Metabolite = c("Propionate", "Glycocholic_acid", "Cholic_acid", "Glycoursodeoxycholic_acid", "Nicotinic_acid", "Indole_acetic_acid"),
  Gene       = c("STAG1", "MSX1", "UCA1", "KMT5A", "SERPINB9", "CFB"), stringsAsFactors = FALSE
)
stopifnot(all(triads$Species %in% colnames(species_log)), all(triads$Metabolite %in% colnames(met_log)),
          all(triads$Gene %in% colnames(host_expr)))


## ---- 1. Edges and node effects ------------------------------------------------------

spearman_edge <- function(a, b) {                                     # a, b named by RNA_sample_id
  ids <- intersect(names(a), names(b))
  ct <- cor.test(a[ids], b[ids], method = "spearman", exact = FALSE)
  c(N = length(ids), Rho = unname(ct$estimate), P = ct$p.value)
}
is_pcr_all <- meta$Response == "pCR"; is_pcr_sp <- sp_meta$Response == "pCR"

edges <- nodes <- list()
for (i in seq_len(nrow(triads))) {
  s <- setNames(species_log[, triads$Species[i]], rownames(species_log))
  m <- setNames(met_log[, triads$Metabolite[i]], rownames(met_log))
  g <- setNames(host_expr[, triads$Gene[i]], rownames(host_expr))
  edges[[i]] <- data.frame(Panel = i, Edge = c("Species-Metabolite", "Metabolite-Gene", "Species-Gene"),
                           rbind(spearman_edge(s, m), spearman_edge(m, g), spearman_edge(s, g)))
  nodes[[i]] <- data.frame(
    Panel = i, Node = c("Species", "Metabolite", "Gene"),
    Label = c(pretty_species(triads$Species[i]), gsub("_", " ", triads$Metabolite[i]), triads$Gene[i]),
    Effect = c(log2((mean(species_raw[is_pcr_sp, triads$Species[i]]) + 1e-6) / (mean(species_raw[!is_pcr_sp, triads$Species[i]]) + 1e-6)),
               log2((median(met_raw[is_pcr_all, triads$Metabolite[i]]) + 1e-6) / (median(met_raw[!is_pcr_all, triads$Metabolite[i]]) + 1e-6)),
               host_deg$LFC[match(triads$Gene[i], host_deg$Gene)])
  )
}
edges <- do.call(rbind, edges); nodes <- do.call(rbind, nodes)
print(cbind(triads, n_edges_P05 = tapply(edges$P < 0.05, edges$Panel, sum)))
save_table(cbind(triads[edges$Panel, ], edges[, -1]), "Fig4E_triad_edges")
save_table(nodes, "Fig4E_triad_nodes")


## ---- 2. Figure 4E ---------------------------------------------------------------------
## Triangle: species bottom-left, metabolite top, gene bottom-right.

tri <- data.frame(Node = c("Species", "Metabolite", "Gene"), x = c(0, 0.72, 1.44), y = c(0, 1.247, 0), label_y = c(-0.3, 1.48, -0.3))
nodes <- merge(nodes, tri, by = "Node")
nodes$Font <- c(Species = "italic", Metabolite = "plain", Gene = "bold")[nodes$Node]
edge_geom <- data.frame(Edge = c("Species-Metabolite", "Metabolite-Gene", "Species-Gene"),
                        x = c(0, 0.72, 0), y = c(0, 1.247, 0), xend = c(0.72, 1.44, 1.44), yend = c(1.247, 0, 0),
                        lx = c(0.28, 1.16, 0.72), ly = c(0.69, 0.69, -0.11), angle = c(60, -60, 0))
edges <- merge(edges, edge_geom, by = "Edge")
edges$Label <- paste0("r=", sprintf("%.2f", edges$Rho), p_stars(edges$P))
edges$Line  <- ifelse(edges$P < 0.05, "solid", "longdash")
edges$Panel <- factor(edges$Panel, levels = seq_len(nrow(triads)), labels = gsub("_", " ", triads$Metabolite))
nodes$Panel <- factor(nodes$Panel, levels = seq_len(nrow(triads)), labels = gsub("_", " ", triads$Metabolite))
effect_lim <- max(1, quantile(abs(nodes$Effect), 0.9))
nodes$Effect_plot <- pmax(-effect_lim, pmin(effect_lim, nodes$Effect))

## the metabolite name is the facet title, so only the species (left, italic) and
## the gene (right, bold) are written under the triangle
corner <- nodes[nodes$Node != "Metabolite", ]
corner$hjust <- ifelse(corner$Node == "Species", 0.7, 0.3)                # species text mostly to the left, gene mostly to the right

p_4e <- ggplot() +
  geom_segment(data = edges, aes(x, y, xend = xend, yend = yend, colour = Rho, linewidth = abs(Rho), linetype = Line), lineend = "round") +
  geom_text(data = edges, aes(lx, ly, label = Label, angle = angle), size = 2.8) +
  geom_point(data = nodes, aes(x, y, fill = Effect_plot), shape = 21, size = 5, colour = "grey30") +
  geom_text(data = corner, aes(x, label_y, label = Label, fontface = Font, hjust = hjust), size = 2.6) +
  scale_colour_gradient2(low = "#3B6EA8", mid = "grey92", high = "#B14A4A", limits = c(-1, 1), name = "Cor.",
                         guide = guide_colourbar(order = 2, direction = "horizontal", title.position = "top",
                                                 barwidth = unit(2.6, "cm"), barheight = unit(0.28, "cm"))) +
  scale_fill_gradient2(low = response_colors[["non-pCR"]], mid = "white", high = response_colors[["pCR"]],
                       limits = c(-effect_lim, effect_lim), name = expression(Log[2] * " FC"),
                       guide = guide_colourbar(order = 1, direction = "horizontal", title.position = "top",
                                               barwidth = unit(2.6, "cm"), barheight = unit(0.28, "cm"))) +
  scale_linewidth(range = c(1, 2.2), guide = "none") +
  scale_linetype_identity() +
  facet_wrap(~Panel, nrow = 1) +
  coord_equal(xlim = c(-0.55, 1.99), ylim = c(-0.5, 1.45), clip = "off") +       # wide x range: room for the corner labels
  theme_void(base_size = 9) +
  theme(strip.text = element_text(size = 9, margin = margin(b = 6)), legend.position = "right", legend.box = "vertical",
        legend.direction = "horizontal", legend.title = element_text(size = 8), legend.text = element_text(size = 7.5),
        legend.spacing.y = unit(3, "mm"), panel.spacing = unit(9, "mm"), plot.margin = margin(4, 4, 6, 30))
save_panel(p_4e, "Fig4E_triads", width = 13, height = 2.5)

message("Done: 17_Fig4E_Species_metabolite_host_triads.R")
