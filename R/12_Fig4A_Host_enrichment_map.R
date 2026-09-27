## =============================================================================
##  12_Fig4A_Host_enrichment_map.R
##  ---------------------------------------------------------------------------
##  PANEL       Figure 4A  enrichment map: response-associated gene sets as
##                         nodes, leading-edge gene overlap as edges, Louvain
##                         communities as labelled hulls
##
##  PURPOSE     Give a bird's-eye view of the whole GSEA result (Hallmark,
##              canonical pathways, GO:BP and colorectal-cancer signatures),
##              grouping redundant gene sets into programmes.
##
##  TWO MODES   refit_gsea <- FALSE (default)  replay from the frozen tables:
##                  the full GSEA family (gsea_global_results_frozen.csv, 8,240
##                  gene sets) is re-screened with the selection rule below, the
##                  edges are recomputed from the leading-edge genes, and the
##                  community membership / labels come from the frozen node
##                  table (Louvain is stochastic).
##              refit_gsea <- TRUE   rerun fgsea on the whole family (needs the
##                  msigdbr package and output/cache/host_deseq2.RData), Louvain
##                  clustering and the keyword labelling rules
##                  (R/_helpers_enrichment_map.R).  Results will differ.
##
##  INPUT       input/host_rnaseq/gsea_global_results_frozen.csv
##              input/host_rnaseq/enrichment_map_nodes_frozen.csv  (160 nodes)
##
##  OUTPUT      output/figures/Fig4A_host_enrichment_map.{svg,png}
##              output/tables/Fig4A_enrichment_map_nodes.csv, Fig4A_enrichment_map_edges.csv
##
##  METHODS     * Family: MSigDB 2026.1.Hs Hallmark + C2 canonical pathways +
##                colorectal C2:CGP signatures + GO:BP, 10-500 genes, fgsea
##                (eps = 0, seed 20260703) on the DESeq2 Wald statistic.
##              * Nodes: FDR < 0.10 and |NES| >= 1.3, or FDR < 0.25 and |NES| >=
##                1.2 for gene sets matching a list of targeted keywords
##                (adenoma, WNT, MYC, G2/M, EMT, hypoxia, TNF/NF-kB ...); at
##                most 80 per direction; gene sets without a biological axis
##                ("Other") excluded.
##              * Edges: leading-edge gene overlap >= 3 genes and Jaccard >= 0.06.
##              * Communities: Louvain (weights = Jaccard); a community is
##                pCR- / non-pCR-dominant when >= 70 % of its nodes and >= 70 % of
##                its evidence weight (|NES| x -log10 FDR) point the same way.
##                Layout: Fruchterman-Reingold (seed 20260703).
##
##  R PACKAGES  dplyr, tidyr, igraph, ggraph, ggforce, ggrepel, stringr
## =============================================================================

source("R/_common.R")
source("R/_helpers_enrichment_map.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(ggforce)
  library(ggrepel)
})
report_versions(c("igraph", "ggraph", "ggforce", "ggrepel"))

refit_gsea <- FALSE      # <-- TRUE = rerun fgsea + Louvain (see header)
seed_map <- 20260703


## ---- 1. The full GSEA family ----------------------------------------------------

if (!refit_gsea) {
  gsea <- read.csv(input_file("host_rnaseq", "gsea_global_results_frozen.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(all(gsea$db_version == "2026.1.Hs"))
  message("Frozen global GSEA table: ", nrow(gsea), " gene sets.")
} else {
  suppressPackageStartupMessages({ library(msigdbr); library(fgsea) })
  load(cache_file("host_deseq2.RData"))                                     # deg, vst (script 10)
  msig <- msigdbr(species = "Homo sapiens")
  gene_col <- intersect(c("gene_symbol", "db_gene_symbol"), names(msig))[1]
  msig <- msig %>% rename(Gene = all_of(gene_col)) %>%
    filter(Gene %in% rownames(vst),
           gs_collection == "H" | (gs_collection == "C2" & grepl("^CP", gs_subcollection)) |
             (gs_collection == "C2" & grepl("CGP", gs_subcollection) & grepl("COLORECTAL|COLON|RECTAL|ADENOMA|\\bCRC\\b", gs_name, ignore.case = TRUE)) |
             (gs_collection == "C5" & gs_subcollection == "GO:BP")) %>%
    distinct(gs_name, Gene, .keep_all = TRUE)
  sizes <- msig %>% count(gs_name, name = "n_genes_total") %>% filter(n_genes_total >= 10, n_genes_total <= 500)
  gene_sets <- split(msig$Gene[msig$gs_name %in% sizes$gs_name], msig$gs_name[msig$gs_name %in% sizes$gs_name])
  ranks <- deg %>% filter(is.finite(stat)) %>% arrange(desc(abs(stat))) %>% distinct(Gene, .keep_all = TRUE)
  ranks <- sort(setNames(ranks$stat, ranks$Gene), decreasing = TRUE)
  set.seed(seed_map)
  gsea <- fgsea(gene_sets, ranks, minSize = 10, maxSize = 500, eps = 0, nproc = 1) %>%
    as.data.frame() %>%
    mutate(gs_name = pathway, display_label = format_pathway_label(gs_name),
           biological_axis = assign_biological_axis(display_label),
           leadingEdge_string = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)),
           Direction = ifelse(NES > 0, "pCR_high", "non_pCR_high"), db_version = "2026.1.Hs") %>%
    select(-leadingEdge) %>% left_join(sizes, by = "gs_name")
}


## ---- 2. Node selection ---------------------------------------------------------------

targeted_regex <- paste(c("adenoma", "colorectal", "colon", "rectal", "\\bCRC\\b", "WNT", "beta.?catenin", "\\bMYC(?:N|L)?\\b",
                          "G2.?M", "E2F", "cell cycle", "mitotic", "DNA replication", "epithelial.*mesenchymal", "\\bEMT\\b",
                          "ECM", "collagen", "oxidative.*stress", "\\bROS\\b", "hypoxia", "TNF", "NF.?KB", "IL.?6", "JAK",
                          "\\bSTAT[0-9AB]*\\b", "\\bTLR[0-9]*\\b", "\\bNOD[0-9]*\\b", "NOD.?LIKE", "cytokine", "chemokine"),
                        collapse = "|")
fdr_cut <- 0.10; nes_cut <- 1.30                    # <-- default node rule
fdr_cut_targeted <- 0.25; nes_cut_targeted <- 1.20  # <-- relaxed rule for targeted keywords
max_per_direction <- 80

nodes <- gsea %>%
  mutate(targeted = grepl(targeted_regex, display_label, ignore.case = TRUE, perl = TRUE)) %>%
  filter(!is.na(padj), !is.na(Direction), biological_axis != "Other",
         (padj < fdr_cut & abs(NES) >= nes_cut) | (targeted & padj < fdr_cut_targeted & abs(NES) >= nes_cut_targeted)) %>%
  group_by(Direction) %>% arrange(padj, desc(abs(NES)), .by_group = TRUE) %>% slice_head(n = max_per_direction) %>% ungroup() %>%
  mutate(node_size = pmin(-log10(padj + 1e-300), 10), evidence_weight = abs(NES) * node_size)
cat("Nodes selected:", nrow(nodes), "(", sum(nodes$Direction == "pCR_high"), "pCR /", sum(nodes$Direction == "non_pCR_high"), "non-pCR )\n")


## ---- 3. Edges: leading-edge gene overlap ----------------------------------------------

leading <- setNames(strsplit(nodes$leadingEdge_string, ";", fixed = TRUE), nodes$gs_name)
pairs <- t(combn(nodes$gs_name, 2))
edges <- data.frame(from = pairs[, 1], to = pairs[, 2])
edges$overlap_n <- mapply(function(a, b) length(intersect(leading[[a]], leading[[b]])), edges$from, edges$to)
edges$jaccard   <- mapply(function(a, b) length(intersect(leading[[a]], leading[[b]])) / length(union(leading[[a]], leading[[b]])), edges$from, edges$to)
edges <- edges[edges$overlap_n >= 3 & edges$jaccard >= 0.06, ]            # <-- edge rule
cat("Edges:", nrow(edges), "\n")


## ---- 4. Communities and labels ---------------------------------------------------------

if (!refit_gsea) {
  frozen_nodes <- read.csv(input_file("host_rnaseq", "enrichment_map_nodes_frozen.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(setequal(nodes$gs_name, frozen_nodes$gs_name))
  i <- match(nodes$gs_name, frozen_nodes$gs_name)
  stopifnot(max(abs(nodes$NES - frozen_nodes$NES[i])) < 1e-10, max(abs(nodes$padj - frozen_nodes$padj[i])) < 1e-10)
  nodes$network_cluster   <- frozen_nodes$network_cluster[i]
  nodes$cluster_label     <- frozen_nodes$cluster_label[i]
  nodes$cluster_direction <- frozen_nodes$cluster_direction[i]
} else {
  set.seed(seed_map)
  g0 <- graph_from_data_frame(edges, vertices = data.frame(name = nodes$gs_name), directed = FALSE)
  nodes$network_cluster <- cluster_louvain(g0, weights = E(g0)$jaccard)$membership[match(nodes$gs_name, V(g0)$name)]
  cl <- nodes %>% group_by(network_cluster) %>%
    summarise(cluster_label = refine_cluster_label(display_label, biological_axis),
              pCR_prop = mean(Direction == "pCR_high"),
              pCR_weighted_prop = sum(evidence_weight[Direction == "pCR_high"]) / sum(evidence_weight), .groups = "drop") %>%
    mutate(cluster_direction = case_when(pCR_prop >= 0.7 & pCR_weighted_prop >= 0.7 ~ "pCR_high",
                                         pCR_prop <= 0.3 & pCR_weighted_prop <= 0.3 ~ "non_pCR_high", TRUE ~ "mixed"))
  nodes <- left_join(nodes, cl, by = "network_cluster")
}

## Community names as printed in the submitted panel (edited from the automatic
## labels; keyed by the automatic label).  Communities not listed here keep
## their automatic label; use NA to leave a community unlabelled.
panel_labels <- c(
  "MYC target\nprogram"                 = "MYC targets",
  "ECM/collagen\nremodeling"            = "Adenoma /\ninflammation and redox programs",
  "Epithelial–mesenchymal\ntransition" = "Epithelial–mesenchymal\ntransition",
  "Translation\ninitiation"             = "Mitochondrial\ntranslation",
  "Cell-cycle\nprogression"             = "Proteasome /\nantigen-processing",
  "MYC-down\n signature"                = "MYC-down\nsignature",
  "MYC-down\nsignature"                 = "MYC-down\nsignature",
  "G2/M checkpoint"                     = "G2/M checkpoint",
  "T-cell / lymphocyte\nactivation"     = "T-cell/lymphocyte\nactivation",
  "Mucus/goblet-cell\nprogram"          = "Mucin glycosylation",
  "Prostaglandin\nresponse"             = NA                    # left unlabelled in the submission
)
nodes$panel_label <- ifelse(nodes$cluster_label %in% names(panel_labels), panel_labels[nodes$cluster_label], nodes$cluster_label)
save_table(nodes %>% select(-any_of("leadingEdge_string")), "Fig4A_enrichment_map_nodes")
save_table(edges, "Fig4A_enrichment_map_edges")


## ---- 5. Figure 4A ------------------------------------------------------------------------

g <- graph_from_data_frame(edges, vertices = nodes %>% rename(name = gs_name) %>% select(name, everything()), directed = FALSE)
set.seed(seed_map)
lay <- create_layout(g, layout = "fr")
lay_df <- as.data.frame(lay)

clusters <- lay_df %>%
  group_by(network_cluster, panel_label, cluster_direction) %>%
  summarise(n_nodes = n(), x = median(x), y = median(y), .groups = "drop") %>%
  filter(n_nodes >= 2)                                            # hulls and labels only for communities of >= 2 gene sets
hull_nodes <- lay_df %>% filter(network_cluster %in% clusters$network_cluster)
labels <- clusters %>% filter(!is.na(panel_label))
x_at <- function(f) min(lay_df$x) + f * diff(range(lay_df$x)); y_at <- function(f) min(lay_df$y) + f * diff(range(lay_df$y))

## Cluster names are placed next to the cluster (offset from its median node
## position, in fractions of the layout width / height) as in the submitted
## panel.  Any cluster not listed here is labelled at its centre.  The seeded
## layout is deterministic, so the offsets stay valid; a refit may need new ones.
label_nudge <- list(
  "G2/M checkpoint"                          = c(-0.11,  0.00),
  "T-cell/lymphocyte\nactivation"            = c(-0.06, -0.01),
  "Mitochondrial\ntranslation"               = c(-0.08,  0.05),
  "Proteasome /\nantigen-processing"         = c( 0.00,  0.08),
  "MYC targets"                              = c( 0.10,  0.00),
  "Mucin glycosylation"                      = c( 0.11,  0.04),
  "MYC-down\nsignature"                      = c( 0.06, -0.08),
  "Epithelial\u2013mesenchymal\ntransition"  = c(-0.11,  0.09),
  "Adenoma /\ninflammation and redox programs" = c(-0.05, -0.07)
)
labels <- labels %>% rowwise() %>%
  mutate(dx = if (panel_label %in% names(label_nudge)) label_nudge[[panel_label]][1] else 0,
         dy = if (panel_label %in% names(label_nudge)) label_nudge[[panel_label]][2] else 0,
         lx = x + dx * diff(range(lay_df$x)), ly = y + dy * diff(range(lay_df$y))) %>% ungroup()

direction_colors <- c(pCR_high = response_colors[["pCR"]], non_pCR_high = response_colors[["non-pCR"]], mixed = "grey70")

p_4a <- ggraph(lay) +
  geom_mark_hull(data = hull_nodes, aes(x, y, group = network_cluster, fill = cluster_direction), inherit.aes = FALSE,
                 alpha = 0.15, colour = "grey55", linewidth = 0.4, concavity = 4, expand = unit(2.5, "mm"), radius = unit(2, "mm"),
                 show.legend = FALSE) +
  geom_edge_link(aes(edge_width = jaccard), colour = "grey70", alpha = 0.45, show.legend = FALSE) +
  scale_edge_width(range = c(0.2, 1.1)) +
  geom_node_point(aes(size = node_size, fill = Direction), shape = 21, colour = "grey25", stroke = 0.3) +
  geom_text(data = labels, aes(lx, ly, label = panel_label), inherit.aes = FALSE, size = 3.3, fontface = "bold", lineheight = 0.85) +
  annotate("text", x = x_at(0.56), y = y_at(0.98), label = "non-pCR-enriched", colour = response_colors[["non-pCR"]], fontface = "bold", size = 4.5) +
  annotate("text", x = x_at(0.47), y = y_at(0.55), label = "pCR-enriched", colour = response_colors[["pCR"]], fontface = "bold", size = 4.5) +
  scale_fill_manual(values = direction_colors, guide = "none") +
  scale_size_continuous(name = expression(-Log[10] ~ FDR), range = c(2, 7.2), limits = c(0, 10), breaks = c(0, 2.5, 5, 7.5, 10)) +
  guides(size = guide_legend(override.aes = list(fill = "grey30"))) +
  scale_x_continuous(expand = expansion(mult = 0.09)) + scale_y_continuous(expand = expansion(mult = 0.07)) +   # room for the labels
  coord_cartesian(clip = "off") +
  theme_void(base_size = 10) +
  theme(legend.position = "inside", legend.position.inside = c(0.97, 0.05), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.4), legend.title = element_text(size = 9))
save_panel(p_4a, "Fig4A_host_enrichment_map", width = 7.2, height = 6.4)

message("Done: 12_Fig4A_Host_enrichment_map.R")
