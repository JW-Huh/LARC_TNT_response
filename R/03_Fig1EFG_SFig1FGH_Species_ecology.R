## =============================================================================
##  03_Fig1EFG_SFig1FGH_Species_ecology.R
##  ---------------------------------------------------------------------------
##  PANELS      Figure 1E  radar chart of six literature-curated microbial groups
##              Figure 1F  dumbbell plot of response-associated species (+ tiles)
##              Figure 1G  Spearman co-abundance network of those species
##              Supplementary Fig. 1F/1G  species heatmaps (baseline / after RT)
##              Supplementary Fig. 1H  RT-induced change of six CRC-associated SGBs
##
##  PURPOSE     Describe which species and functional guilds differ between
##              pCR and non-pCR at baseline, how those species co-vary, and how
##              prespecified CRC-associated taxa change after radiotherapy.
##
##  INPUT       input/microbiome/metaphlan_species.csv   species abundance (%)
##              input/microbiome/metaphlan_sgb.csv       SGB abundance (%) (S1H)
##              input/microbiome/microbial_groups.csv    curated group memberships
##                                                        (Supplementary Table 1)
##              input/metadata/stool_samples.csv, patients.csv
##
##  OUTPUT      output/figures/Fig1E_microbial_group_radar, Fig1F_species_dumbbell,
##              Fig1G_species_network, SFig1F_species_heatmap_baseline,
##              SFig1G_species_heatmap_afterRT, SFig1H_CRC_SGB_change  (.svg/.png)
##              output/figures/Fig1G_species_network.graphml   (Cytoscape input)
##              output/tables/Fig1E_group_scores.csv, Fig1F_species_statistics.csv,
##              Fig1G_network_edges.csv, SFig1H_paired_changes.csv
##
##  METHODS
##   * Group scores (1E): for each sample the relative abundances of the member
##     species are summed, log10(sum + 1) transformed and divided by the group's
##     90th percentile across the 26 samples (x100).  The radar shows the mean
##     score of each response group.  Memberships follow the figure-producing
##     code: the ICB and CRC groups are the union of the Table S1A list and the
##     Table S1B annotations (22 and 21 taxa).
##   * Species screen (1F, 1G, S1F): prevalence >= 20 % of baseline samples,
##     two-sided Wilcoxon rank-sum P < 0.1, and a resolved species name (SGBs
##     named "GGB..." are excluded).  Points = log10(group mean abundance + 1e-4);
##     log2FC = log2((mean_pCR + 1e-4) / (mean_nonpCR + 1e-4)).
##   * Network (1G): Spearman correlations between all pairs of screened
##     species (26 baseline samples); edges with nominal P < 0.01 are kept.
##     Node positions reproduce the submitted Cytoscape layout (see below).
##   * Heatmaps (S1F/S1G): log10(abundance + pseudocount), pseudocount = half
##     of the smallest positive value; row z-scores; rows clustered with
##     Manhattan distance / Ward.D2, samples ordered within each response group
##     by Euclidean / complete clustering.  After-RT rows are the baseline-
##     selected species with prevalence > 20 % after RT.
##   * S1H: for the 16 patients with paired samples, log10((after + 1e-5) /
##     (before + 1e-5)) per SGB; mean and normal-approximation 95 % CI.
##
##  R PACKAGES  dplyr, tidyr, ggplot2, ggnewscale, patchwork, ComplexHeatmap,
##              circlize, igraph, ggraph
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggnewscale)
  library(ComplexHeatmap)
  library(circlize)
  library(igraph)
  library(ggraph)
})
report_versions(c("dplyr", "ggplot2", "ComplexHeatmap", "igraph", "ggraph"))
set_heatmap_options()                                   # legend font sizes (R/_common.R)


## ---- 1. Data --------------------------------------------------------------------

samples  <- read_stool_samples()                                   # 42 samples
baseline <- samples[samples$Timepoint == "Before", ]                # 26
species_all <- read_metaphlan("species", samples = samples$SampleID, drop_empty = FALSE)
sp_base <- species_all[baseline$SampleID, ]                         # 26 x 1227 (%)

groups <- read.csv(input_file("microbiome", "microbial_groups.csv"), stringsAsFactors = FALSE)
member <- function(list_name) unique(groups$Species[groups$List %in% list_name])
stopifnot(all(groups$Species %in% colnames(species_all)))

## Six radar groups (Figure 1E).  ICB and CRC unite the S1A list with the S1B
## annotations, exactly as in the code that produced the submitted panel.
radar_groups <- list(
  "Cardiometabolic health\n-associated bacteria" = member("health_S1A"),
  "ICB-favorable\nbacteria"    = member(c("icb_S1A", "icb_S1B")),
  "Sulfidogenic\nbacteria"     = member("sulfidogenic_S1A"),
  "Oral-derived\nCRC signatures" = member(c("crc_S1A", "crc_S1B")),
  "Mucin\ndegraders"           = member("mucin_S1A"),
  "Butyrate-\nproducing guild" = member(c("butyrate_S1A", "butyrate_producer_S1B", "butyrate_crossfeeder_S1B"))
)
cat("Radar group sizes:", paste(sapply(radar_groups, length), collapse = "/"), "\n")

## Tile annotations (Figure 1F/1G): Table S1B categories; a species that is both
## a direct producer and a cross-feeder is shown as a direct producer.
tile_producer <- member("butyrate_producer_S1B")
tile_feeder   <- member("butyrate_crossfeeder_S1B")
tile_icb      <- member("icb_S1B")
tile_crc      <- member("crc_S1B")
oral_crc      <- member("crc_S1A")                                  # white circle in the CRC tile


## ---- 2. Figure 1E: radar chart of group scores --------------------------------------

score_of <- function(sp) {
  v <- log10(rowSums(sp_base[, sp, drop = FALSE]) + 1)
  100 * v / quantile(v, 0.90)                                       # <-- scaling quantile
}
scores <- data.frame(SampleID = baseline$SampleID, Response = baseline$Response,
                     sapply(radar_groups, score_of), check.names = FALSE)
radar_mean <- scores %>%
  pivot_longer(-c(SampleID, Response), names_to = "Group", values_to = "Score") %>%
  group_by(Response, Group) %>% summarise(Score = mean(Score), .groups = "drop") %>%
  mutate(Group = factor(Group, levels = names(radar_groups)))
save_table(radar_mean, "Fig1E_group_scores")

## polar coordinates: first group at 12 o'clock, then clockwise
angles <- setNames(seq(pi / 2, pi / 2 - 2 * pi, length.out = 7)[1:6], names(radar_groups))
radar_mean <- radar_mean %>%
  mutate(angle = angles[as.character(Group)], x = Score * cos(angle), y = Score * sin(angle)) %>%
  arrange(Response, angle)
radar_closed <- bind_rows(radar_mean, radar_mean %>% group_by(Response) %>% slice(1) %>% ungroup())
grid_levels <- c(20, 40, 60, 80, 100)
radar_grid <- do.call(rbind, lapply(grid_levels, function(r)
  data.frame(r = r, angle = c(angles, angles[1]), x = r * cos(c(angles, angles[1])), y = r * sin(c(angles, angles[1])))))
radar_axis <- data.frame(xend = 100 * cos(angles), yend = 100 * sin(angles))
radar_label <- data.frame(label = names(angles), x = 122 * cos(angles), y = 118 * sin(angles))

p_radar <- ggplot() +
  geom_path(data = radar_grid, aes(x, y, group = r), colour = "grey85", linewidth = 0.4) +
  geom_segment(data = radar_axis, aes(x = 0, y = 0, xend = xend, yend = yend), colour = "grey85", linewidth = 0.4) +
  geom_polygon(data = radar_closed, aes(x, y, group = Response, fill = Response), alpha = 0.4) +
  geom_path(data = radar_closed, aes(x, y, group = Response, colour = Response), linewidth = 1) +
  geom_point(data = radar_mean, aes(x, y, colour = Response), size = 2) +
  geom_text(data = radar_label, aes(x, y, label = label), size = 3.6, lineheight = 0.9) +
  scale_fill_manual(values = response_colors, name = "TRG") +
  scale_colour_manual(values = response_colors, name = "TRG") +
  coord_equal(clip = "off") +
  theme_void() +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.95), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
        plot.margin = margin(15, 40, 15, 40))
save_panel(p_radar, "Fig1E_microbial_group_radar", width = 5.6, height = 5)


## ---- 3. Species screen and Figure 1F dumbbell plot ---------------------------------

long <- data.frame(SampleID = rownames(sp_base), sp_base, check.names = FALSE) %>%
  pivot_longer(-SampleID, names_to = "Species", values_to = "Abundance") %>%
  left_join(baseline[, c("SampleID", "Response")], by = "SampleID")

species_stats <- long %>%
  group_by(Species) %>%
  summarise(
    prevalence = mean(Abundance > 0),
    mean_pCR = mean(Abundance[Response == "pCR"]),
    mean_nonpCR = mean(Abundance[Response == "non-pCR"]),
    P = wilcox_p(Abundance, Response),
    .groups = "drop"
  ) %>%
  mutate(
    log2FC = log2((mean_pCR + 1e-4) / (mean_nonpCR + 1e-4)),
    Butyrate = case_when(Species %in% tile_producer ~ "Direct butyrate producer",
                         Species %in% tile_feeder ~ "Butyrate cross-feeder", TRUE ~ "Unassigned"),
    ICB = ifelse(Species %in% tile_icb, "ICB-favorable bacteria", "Unassigned"),
    CRC = ifelse(Species %in% tile_crc, "CRC signatures", "Unassigned"),
    Oral = Species %in% oral_crc
  )

## Screening thresholds  <-- change here to alter which species are displayed
selected <- species_stats %>%
  filter(prevalence >= 0.20, P < 0.10, !grepl("GGB", Species)) %>%
  arrange(desc(log10(mean_pCR + 1e-4)))
cat("Screened species:", nrow(selected), "\n")
save_table(selected, "Fig1F_species_statistics")

## taxon_plotmath() (R/_common.R) writes the binomial in italics and strain /
## SGB codes upright, e.g. italic("Ruminococcus")~"sp. AF13-28"
species_order <- rev(selected$Species)                              # highest pCR mean at the top
sel <- selected %>% mutate(
  Species = factor(Species, levels = species_order),
  label = taxon_plotmath(as.character(Species)),
  stars = p_stars(P, cut = c(0.01, 0.05), symbols = c("**", "*"))
)

## left block: species names
p_names <- ggplot(sel, aes(x = 1, y = Species)) +
  geom_text(aes(label = label), parse = TRUE, hjust = 1, size = 3.1) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  theme_void()

## middle block: category tiles + P value + log2FC tiles
tiles_cat <- sel %>% select(Species, Butyrate, ICB, CRC) %>%
  pivot_longer(-Species, names_to = "Column", values_to = "Group") %>%
  mutate(Column = factor(Column, levels = c("Butyrate", "ICB", "CRC")))
p_tiles <- ggplot() +
  geom_tile(data = tiles_cat, aes(Column, Species, fill = Group), colour = "grey70", linewidth = 0.3) +
  scale_fill_manual(values = group_tile_colors, name = "Group", breaks = names(group_tile_colors)[1:4]) +
  geom_point(data = filter(sel, Oral), aes(x = "CRC", y = Species), shape = 21, fill = "white", colour = "white", size = 2.2) +
  new_scale_fill() +
  geom_tile(data = sel, aes(x = "-log10(p)", y = Species, fill = -log10(P)), colour = "grey70", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(1, 2.8), oob = scales::squish,
                      name = expression(-Log[10](italic(p)))) +
  geom_text(data = sel, aes(x = "-log10(p)", y = Species, label = stars), colour = "white", size = 3.5, vjust = 0.75) +
  new_scale_fill() +
  geom_tile(data = sel, aes(x = "Log2FC", y = Species, fill = log2FC), colour = "grey70", linewidth = 0.3) +
  scale_fill_gradient2(low = response_colors[["non-pCR"]], mid = "white", high = response_colors[["pCR"]],
                       limits = c(-5, 5), oob = scales::squish, name = expression(Log[2]*FC)) +
  scale_x_discrete(limits = c("Butyrate", "ICB", "CRC", "-log10(p)", "Log2FC"), expand = c(0, 0),
                   labels = c("Butyrate", "ICB", "CRC", expression(italic(p)*"-value"), "LFC")) +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))

## right block: dumbbells
dumb <- sel %>% mutate(pCR = log10(mean_pCR + 1e-4), `non-pCR` = log10(mean_nonpCR + 1e-4))
p_dumb <- ggplot(dumb) +
  geom_segment(aes(y = Species, yend = Species, x = `non-pCR`, xend = pCR), colour = "grey70", linewidth = 0.8) +
  geom_point(aes(x = `non-pCR`, y = Species, colour = "non-pCR"), size = 3) +
  geom_point(aes(x = pCR, y = Species, colour = "pCR"), size = 3) +
  scale_colour_manual(values = response_colors, breaks = response_levels, name = "TRG") +
  scale_y_discrete(expand = expansion(add = 0.6)) +
  labs(x = expression(Log[10]~"(mean abundance + 1e"^-4*")"), y = NULL) +
  theme_tnt(base_size = 10, legend = "right") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_1f <- (p_names | p_tiles | p_dumb) + plot_layout(widths = c(2.6, 1.7, 4.2), guides = "collect") &
  theme(legend.position = "right", legend.justification = "top", legend.title = element_text(size = 10),
        legend.text = element_text(size = 8.5), legend.key.size = unit(0.45, "cm"))
save_panel(p_1f, "Fig1F_species_dumbbell", width = 9.5, height = 6.5)


## ---- 4. Figure 1G: Spearman co-abundance network -----------------------------------

net_species <- selected$Species
pairs <- t(combn(net_species, 2))
edges <- data.frame(from = pairs[, 1], to = pairs[, 2], rho = NA_real_, P = NA_real_)
for (i in seq_len(nrow(edges))) {
  ct <- suppressWarnings(cor.test(sp_base[, edges$from[i]], sp_base[, edges$to[i]], method = "spearman", exact = FALSE))
  edges$rho[i] <- unname(ct$estimate); edges$P[i] <- ct$p.value
}
edges$FDR <- p.adjust(edges$P, "BH")
edges <- edges[edges$P < 0.01, ]                                    # <-- edge threshold
edges$sign <- ifelse(edges$rho > 0, "Positive", "Negative")
save_table(edges, "Fig1G_network_edges")

nodes <- data.frame(name = union(edges$from, edges$to))
nodes$mean_abundance <- colMeans(sp_base[, nodes$name, drop = FALSE])
nodes$Group <- with(species_stats[match(nodes$name, species_stats$Species), ],
                    ifelse(Butyrate != "Unassigned", Butyrate, ifelse(CRC != "Unassigned", "CRC signatures", "Unassigned")))
nodes$label <- taxon_plotmath(nodes$name)                          # italic binomial (plotmath)
cat("Network:", nrow(nodes), "nodes,", nrow(edges), "edges (", sum(edges$rho < 0), "negative )\n")

g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
write_graph(g, file.path(paths$figures, "Fig1G_species_network.graphml"), format = "graphml")

## Node coordinates of the submitted figure (manually arranged in Cytoscape and
## read back from the vector PDF; in points, y increases downwards).  Any node
## absent from this table gets a Fruchterman-Reingold position instead.
## Set use_submitted_layout <- FALSE to use a purely algorithmic layout.
use_submitted_layout <- TRUE
submitted_xy <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
Species                                 x       y
Lachnospiraceae_bacterium_CLA_AA_H244   214.8   449.1
Bifidobacterium_dentium                 161.0   454.4
Schaalia_SGB17154                       228.7   476.3
Veillonella_dispar                      149.2   487.0
Rothia_mucilaginosa                     274.3   490.5
Lancefieldella_parvula                  233.9   506.4
Clostridiales_bacterium_KLE1615         135.9   508.2
Allisonella_histaminiformans            106.5   518.3
Propionibacterium_acidifaciens          258.5   528.6
Ruminococcus_callidus                    87.1   529.9
Actinomyces_graevenitzii                207.9   538.2
Roseburia_inulinivorans                 112.7   543.7
Ruminococcus_sp_AF13_28                 150.6   552.0
Solobacterium_SGB6828                   238.1   556.7
Clostridium_SGB4750                      74.8   557.9
Eisenbergiella_tayi                     137.7   566.7
Ruminococcus_sp_AF41_9                  105.9   573.9
Turicibacter_sanguinis                  244.0   585.5
Eubacterium_rectale                     122.0   594.8
Lacrimispora_amygdalina                  73.1   606.3
Clostridium_saudiense                   280.2   610.0
Gemella_sanguinis                       223.5   616.2")

set.seed(123)
xy <- layout_with_fr(g, weights = abs(E(g)$rho))
if (use_submitted_layout) {
  hit <- match(V(g)$name, submitted_xy$Species)
  if (all(!is.na(hit))) {
    xy <- cbind(submitted_xy$x[hit], -submitted_xy$y[hit])
  } else {
    warning("Node set differs from the submitted figure; using the algorithmic layout.")
  }
}
lay <- create_layout(g, layout = "manual", x = xy[, 1], y = xy[, 2])

p_net <- ggraph(lay) +
  geom_edge_link(aes(colour = sign), width = 1.1, alpha = 0.9) +
  scale_edge_colour_manual(values = c(Positive = "#7FB3D5", Negative = "#E08A83"), name = "Correlation") +
  geom_node_point(aes(fill = Group, size = mean_abundance), shape = 21, colour = "white", stroke = 0.4) +
  ## node size = mean baseline abundance (no legend, as submitted); labels are
  ## pushed off the nodes with ggrepel (seed fixed for a reproducible placement)
  geom_node_text(aes(label = label), parse = TRUE, repel = TRUE, seed = 1, size = 2.9, colour = "black",
                 point.padding = unit(0.3, "lines"), box.padding = unit(0.15, "lines"), segment.colour = NA) +
  scale_fill_manual(values = group_tile_colors, name = "Group", breaks = names(group_tile_colors)[1:4]) +
  scale_size_continuous(range = c(5, 15), guide = "none") +
  coord_equal(clip = "off") +
  theme_void() +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8),
        legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.3),
        plot.margin = margin(10, 30, 10, 10))
save_panel(p_net, "Fig1G_species_network", width = 9, height = 5.5)


## ---- 5. Supplementary Figure 1F / 1G: species heatmaps -----------------------------


species_heatmap <- function(timepoint, name) {
  meta <- samples[samples$Timepoint == timepoint, ]
  x <- t(species_all[meta$SampleID, net_species, drop = FALSE])
  if (timepoint == "Ongoing") x <- x[rowMeans(x > 0) > 0.20, , drop = FALSE]     # after-RT prevalence > 20 %
  pseudo <- min(x[x > 0]) / 2
  z <- t(scale(t(log10(x + pseudo))))
  p <- apply(x, 1, function(v) wilcox_p(v, meta$Response, exact = FALSE))

  ## samples ordered within response group (pCR first) by Euclidean / complete clustering
  col_order <- order_samples_within_response(z, meta)
  meta <- meta[match(col_order, meta$SampleID), ]
  z <- z[, col_order]

  top <- response_top_annotation(meta)                 # TRG score + TRG tiles (R/_common.R)
  right <- rowAnnotation(
    `p-value` = anno_simple(-log10(p), col = colorRamp2(c(0, 3), c("white", "#2166AC")),
                            pch = p_stars(p, c(0.05, 0.10), c("**", "*")), pt_gp = gpar(col = "white", fontsize = 9),
                            gp = gpar(col = "grey40", lwd = 0.5), width = unit(4, "mm")),
    annotation_name_gp = gpar(fontsize = 9, fontface = "italic"), annotation_name_rot = 90
  )
  ht <- Heatmap(
    z, name = "Z-score", col = colorRamp2(c(-2, 0, 2), zscore_colors),
    cluster_columns = FALSE, column_split = meta$Response, column_title = NULL, column_gap = unit(1.5, "mm"),
    clustering_distance_rows = "manhattan", clustering_method_rows = "ward.D2",
    row_dend_side = "left", row_dend_width = unit(1.2, "cm"),
    show_column_names = FALSE, row_names_side = "right",
    row_labels = pretty_species(rownames(z)), row_names_gp = gpar(fontsize = 8.5, fontface = "italic"),
    rect_gp = gpar(col = "grey40", lwd = 0.5), top_annotation = top, right_annotation = right,
    heatmap_legend_param = list(direction = "horizontal", at = -2:2, title_position = "topcenter")
  )
  lgd_p <- Legend(title = expression(-log[10](italic(p))), col_fun = colorRamp2(c(0, 3), c("white", "#2166AC")),
                  at = 0:3, direction = "horizontal", title_position = "topcenter")
  dev <- open_panel(name, width = ifelse(timepoint == "Before", 7.8, 6.4), height = ifelse(timepoint == "Before", 8, 6.2))
  draw_panel(dev, draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                       merge_legends = TRUE, annotation_legend_list = list(lgd_p)))
  invisible(data.frame(Species = rownames(z), P = p))
}
species_heatmap("Before",  "SFig1F_species_heatmap_baseline")
species_heatmap("Ongoing", "SFig1G_species_heatmap_afterRT")


## ---- 6. Supplementary Figure 1H: paired change of CRC-associated SGBs -------------

crc_sgbs <- c("Hungatella_hathewayi|SGB4742", "Fusobacterium_nucleatum|SGB6011", "Gemella_morbillorum|SGB7295",
              "Solobacterium_moorei|SGB6826", "Peptostreptococcus_stomatis|SGB748", "Parvimonas_micra|SGB6653")
sgb_all <- read_metaphlan("sgb", samples = samples$SampleID, drop_empty = FALSE)
stopifnot(all(crc_sgbs %in% colnames(sgb_all)))

paired <- merge(samples[samples$Timepoint == "Before", c("SubjectID", "SampleID")],
                samples[samples$Timepoint == "Ongoing", c("SubjectID", "SampleID")],
                by = "SubjectID", suffixes = c("_before", "_after"))
cat("Patients with paired metagenomes:", nrow(paired), "\n")
fc <- log10((sgb_all[paired$SampleID_after, crc_sgbs] + 1e-5) / (sgb_all[paired$SampleID_before, crc_sgbs] + 1e-5))
change <- data.frame(
  SGB = crc_sgbs, mean = colMeans(fc), se = apply(fc, 2, sd) / sqrt(nrow(fc))
) %>% mutate(lower = mean - 1.96 * se, upper = mean + 1.96 * se,
             label = sub("\\|", " | ", pretty_species(SGB))) %>%
  arrange(desc(mean))
change$label <- factor(change$label, levels = rev(change$label))
save_table(change, "SFig1H_paired_changes")

p_1h <- ggplot(change, aes(mean, label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.25, linewidth = 0.5, orientation = "y") +
  geom_point(shape = 21, size = 3.2, fill = "#F5BE1E", colour = "black") +
  annotate("text", x = -0.05, y = nrow(change) + 0.9, hjust = 1, size = 3, label = "Decreased after RT") +
  annotate("text", x = 0.05, y = nrow(change) + 0.9, hjust = 0, size = 3, label = "Increased") +
  annotate("segment", x = -0.1, xend = -0.95, y = nrow(change) + 0.9, yend = nrow(change) + 0.9,
           arrow = arrow(length = unit(0.15, "cm"))) +
  annotate("segment", x = 0.65, xend = 0.95, y = nrow(change) + 0.9, yend = nrow(change) + 0.9,
           arrow = arrow(length = unit(0.15, "cm"))) +
  coord_cartesian(ylim = c(0.5, nrow(change) + 1.2), clip = "off") +
  labs(x = expression(Log[10]~"Fold Change (After RT / Baseline)"), y = NULL) +
  theme_tnt(base_size = 10, legend = "none") +
  theme(axis.text.y = element_text(face = "italic"))
save_panel(p_1h, "SFig1H_CRC_SGB_change", width = 5.6, height = 2.8)

message("Done: 03_Fig1EFG_SFig1FGH_Species_ecology.R")
