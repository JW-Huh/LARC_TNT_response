## =============================================================================
##  13_Fig4C_Species_host_correlations.R
##  ---------------------------------------------------------------------------
##  PANEL       Figure 4C  Spearman correlations between baseline gut species and
##                         18 response-associated tumour genes (20 patients with
##                         both a baseline metagenome and a baseline biopsy)
##
##  PURPOSE     Ask whether the abundance of individual gut species tracks the
##              expression of host genes drawn from the enriched programmes of
##              Figure 4B (adenoma/EMT, gut homing, T-cell differentiation,
##              lymphocyte immunity / TNF signalling).
##
##  INPUT       output/cache/host_deseq2.RData        vst, deg, rna_meta (script 10)
##              output/cache/host_enrichment.RData    gsea_annotation (script 11)
##              input/microbiome/metaphlan_species.csv
##
##  OUTPUT      output/figures/Fig4C_species_host_correlations.{svg,png}
##              output/tables/Fig4C_species_host_correlations.csv  (all tested pairs)
##              output/tables/Fig4C_displayed_species.csv
##
##  METHODS
##   * Host panel: 18 prespecified genes in four blocks (see `host_genes`).
##     The top annotation marks, for up to 10 enriched gene sets (Figure 4B
##     catalogue, FDR < 0.1, greedy selection of sets that add >= 30 % new
##     leading-edge genes), whether a gene is in the set's leading edge.
##   * Species screen (20 matched patients): prevalence >= 30 % and >= 6
##     non-zero samples; "unclear" names (GGB/SGB/sp./bacterium ...) are dropped
##     unless the species is on the a-priori priority list (CRC-, oral- or
##     immune-relevant taxa); Weissella confusa and Anaerostipes hadrus are
##     excluded a priori; species whose rank profile is practically identical
##     (Spearman >= 0.995) to an already kept species are collapsed.
##   * Correlation: Spearman rho between log10(relative abundance + 1e-6) and
##     VST expression, normal-approximation P, n >= 8.
##   * Display: species with >= 1 gene at P < 0.05 or >= 2 genes at P < 0.10,
##     ranked by number of significant genes (E. rectale and R. inulinivorans
##     always shown; B. fragilis not shown), capped at 31 rows.  Rows: Euclidean
##     clustering of the rho profiles; columns clustered within each host block.
##     This is a curated display rule, not an FDR-controlled screen (see the
##     BH FDR column in the results table).
##
##  R PACKAGES  dplyr, tidyr, stringr, ComplexHeatmap, circlize
##
##  NOTE        When run with Rscript, ComplexHeatmap measures label widths on a
##              temporary pdf device, which cannot encode the en dash / kappa of
##              "PI3K–AKT signaling with NF-κB" and prints two harmless
##              'conversion failure ... mbcsToSbcs' warnings.  The SVG and PNG
##              (svglite / cairo) contain the correct characters.
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
})
report_versions(c("ComplexHeatmap", "circlize", "dplyr"))
set_heatmap_options()

load(cache_file("host_deseq2.RData"))                      # vst, deg, rna_meta
load(cache_file("host_enrichment.RData"))                  # gsea_annotation
stopifnot(nrow(rna_meta) == 20, identical(colnames(vst), rna_meta$RNA_sample_id))
is_pcr <- rna_meta$Response == "pCR"


## ---- 1. Species table of the 20 matched patients ------------------------------------

species_tab <- read.csv(input_file("microbiome", "metaphlan_species.csv"), check.names = FALSE, stringsAsFactors = FALSE)
taxonomy <- data.frame(
  Species = sub(".*s__", "", species_tab$clade_name),
  Phylum  = str_match(species_tab$clade_name, "p__([^|]+)")[, 2],
  Family  = str_match(species_tab$clade_name, "f__([^|]+)")[, 2], stringsAsFactors = FALSE
)
species <- read_metaphlan("species", samples = rna_meta$SampleID, drop_empty = TRUE)   # 20 x species (%)
rownames(species) <- rna_meta$RNA_sample_id


## ---- 2. Host gene panel and gene-set membership annotation ----------------------------

host_genes <- c(
  "MSX1", "NOTUM", "UCA1", "NKD1", "HTRA1", "PLOD1",        # Adenoma signature + EMT
  "NINJ1", "CX3CL1", "ADAM8", "ICAM1",                        # Gut homing + leukocyte trafficking
  "LIG4", "CLPTM1", "METTL14", "ZBTB7B",                      # T-cell differentiation
  "PIK3CB", "CRK", "KLRK1", "RAB27A"                          # Lymphocyte immunity + TNF signalling
)
host_block <- factor(rep(c("Adenoma signature + EMT", "Gut homing + leukocyte trafficking",
                           "T-cell differentiation", "Lymphocyte immunity + TNF signaling"), c(6, 4, 4, 4)),
                     levels = c("Adenoma signature + EMT", "Gut homing + leukocyte trafficking",
                                "T-cell differentiation", "Lymphocyte immunity + TNF signaling"))
names(host_block) <- host_genes
stopifnot(all(host_genes %in% rownames(vst)))

## Greedy choice of <= 10 enriched gene sets whose leading edges cover distinct genes
gs <- gsea_annotation %>%
  filter(pass_display_filter, padj < 0.10) %>%
  mutate(edge = str_split(coalesce(leadingEdge, ""), ";"),
         n_in_vst = vapply(edge, function(e) sum(e %in% rownames(vst)), integer(1))) %>%
  filter(n_in_vst >= 5) %>% arrange(padj, desc(abs(NES)), desc(n_in_vst))
keep <- logical(nrow(gs)); covered <- character(0)
for (i in seq_len(nrow(gs))) {
  if (sum(keep) >= 10) break
  new_frac <- length(setdiff(gs$edge[[i]], covered)) / max(1, length(gs$edge[[i]]))
  if (sum(keep) == 0 || new_frac >= 0.30 || sum(keep) < 6) { keep[i] <- TRUE; covered <- union(covered, gs$edge[[i]]) }
}
gs <- gs[keep, ]

## Annotation rows = the selected gene sets, named as in Figure 4B (display_label
## was harmonised in script 11); one shorter name as printed in the submitted panel
short_label <- c("Colorectal adenoma-associated signature" = "Colorectal adenoma signature")
membership <- gs %>% select(display_label, padj, edge) %>% unnest(edge) %>%
  filter(edge %in% host_genes) %>% distinct(display_label, padj, Gene = edge) %>%
  mutate(label = ifelse(display_label %in% names(short_label), short_label[display_label], display_label)) %>%
  arrange(padj)
category_mat <- matrix("Not included", length(unique(membership$label)), length(host_genes),
                       dimnames = list(unique(membership$label), host_genes))
for (k in seq_len(nrow(membership))) category_mat[membership$label[k], membership$Gene[k]] <- "Included"
category_mat <- category_mat[!duplicated(apply(category_mat, 1, paste, collapse = "")), , drop = FALSE]   # identical rows collapsed

host_stat <- deg %>% filter(Gene %in% host_genes) %>% distinct(Gene, .keep_all = TRUE) %>%
  transmute(Gene, LFC = log2FoldChange, P = P, neglog10P = pmin(-log10(pmax(P, .Machine$double.xmin)), 10), stars = p_stars(P))
host_stat <- host_stat[match(host_genes, host_stat$Gene), ]


## ---- 3. Species screen ---------------------------------------------------------------

priority_regex <- regex(paste(c(
  "Eubacterium[_ ]rectale", "Roseburia[_ ]inulinivorans", "Faecalibacterium[_ ]prausnitzii", "Coprococcus[_ ]catus",
  "Clostridium[_ ]bolteae", "^Fusobacterium([_ ]|$)", "Parvimonas[_ ]micra", "Hungatella[_ ]hathewayi", "Gemella[_ ]morbillorum",
  "Peptostreptococcus[_ ]stomatis", "(Segatella|Prevotella)[_ ]copri", "Bacteroides[_ ]fragilis", "Clostridium[_ ]symbiosum",
  "Escherichia[_ ]coli", "^Campylobacter([_ ]|$)", "Ruminococcus[_ ]torques", "(Ruminococcus|Mediterraneibacter)[_ ]gnavus",
  "Prevotella[_ ]stercorea", "Bacteroides[_ ]thetaiotaomicron", "Weissella[_ ]cibaria", "Fusicatenibacter[_ ]saccharivorans",
  "^Porphyromonas([_ ]|$)", "^Bifidobacterium([_ ]|$)", "Streptococcus[_ ]salivarius", "Streptococcus[_ ]sanguinis", "^Veillonella([_ ]|$)"),
  collapse = "|"), ignore_case = TRUE)                         # <-- a-priori priority taxa
excluded_regex <- regex("Weissella[_ ]confusa|Anaerostipes[_ ]hadrus", ignore_case = TRUE)
unclear_regex  <- regex("GGB|SGB|CAG|UBA|MAG|uncultured|metagenome|_sp_|bacterium|oral_taxon", ignore_case = TRUE)
forced_regex   <- regex("Eubacterium[_ ]rectale|Roseburia[_ ]inulinivorans", ignore_case = TRUE)   # always displayed

sp_stat <- data.frame(Species = colnames(species),
                      mean_pCR = colMeans(species[is_pcr, ]), mean_nonpCR = colMeans(species[!is_pcr, ]),
                      prevalence = colMeans(species > 0), nonzero_n = colSums(species > 0), mean_abund = colMeans(species),
                      P = apply(species, 2, function(v) wilcox_p(v, rna_meta$Response, exact = FALSE)), row.names = NULL) %>%
  mutate(logFC = log2((mean_pCR + 1e-6) / (mean_nonpCR + 1e-6)),
         priority = str_detect(Species, priority_regex), excluded = str_detect(Species, excluded_regex),
         unclear = str_detect(Species, unclear_regex)) %>%
  filter(!excluded, prevalence >= 0.30, nonzero_n >= 6, !unclear | priority) %>%          # <-- species filters
  arrange(desc(priority), P, desc(prevalence), desc(mean_abund), desc(abs(logFC)))

## collapse species with (near-)identical rank profiles; priority species are never collapsed
rank_mat <- apply(species[, sp_stat$Species], 2, rank)
kept <- character(0)
for (sp in sp_stat$Species) {
  if (sd(rank_mat[, sp]) == 0) next
  if (sp_stat$priority[sp_stat$Species == sp] || length(kept) == 0) { kept <- c(kept, sp); next }
  sim <- suppressWarnings(cor(rank_mat[, sp], rank_mat[, kept, drop = FALSE], method = "spearman")); sim[!is.finite(sim)] <- 0
  if (max(sim) < 0.995) kept <- c(kept, sp)
}
sp_stat <- sp_stat %>% filter(Species %in% kept) %>% slice_head(n = 200)
cat("Species entering the correlation screen:", nrow(sp_stat), "\n")


## ---- 4. Spearman correlations species x gene ---------------------------------------------

x_sp <- log10(species[, sp_stat$Species] + 1e-6)
y_ge <- t(vst[host_genes, rna_meta$RNA_sample_id])
cors <- expand_grid(Species = sp_stat$Species, Gene = host_genes)
cors[, c("rho", "P")] <- t(mapply(function(s, g) {
  x <- x_sp[, s]; y <- y_ge[, g]
  if (sum(is.finite(x) & is.finite(y)) < 8 || sd(x) == 0 || sd(y) == 0) return(c(NA, NA))
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)); c(unname(ct$estimate), ct$p.value)
}, cors$Species, cors$Gene))
cors$FDR <- p.adjust(cors$P, "BH")
save_table(cors, "Fig4C_species_host_correlations")


## ---- 5. Species to display ---------------------------------------------------------------

sp_rank <- cors %>% group_by(Species) %>%
  summarise(best_p = min(c(P, 1), na.rm = TRUE), n_sig_genes = n_distinct(Gene[!is.na(P) & P < 0.05]),
            n_p10_genes = n_distinct(Gene[!is.na(P) & P < 0.10]), .groups = "drop") %>%
  left_join(cors %>% mutate(Block = host_block[Gene]) %>% group_by(Species, Block) %>%
              summarise(block_sig = sum(P < 0.05, na.rm = TRUE), block_p10 = sum(P < 0.10, na.rm = TRUE),
                        block_consistency = max(mean(rho >= 0, na.rm = TRUE), mean(rho <= 0, na.rm = TRUE)),
                        block_abs_rho = abs(mean(rho, na.rm = TRUE)), .groups = "drop") %>%
              group_by(Species) %>% arrange(desc(block_sig), desc(block_p10), desc(block_consistency), desc(block_abs_rho), .by_group = TRUE) %>%
              slice_head(n = 1) %>% ungroup() %>% select(-Block), by = "Species") %>%
  left_join(sp_stat, by = "Species") %>%
  mutate(forced = str_detect(Species, forced_regex),
         tier = case_when(forced ~ 0L, n_sig_genes >= 2 ~ 1L, n_sig_genes >= 1 ~ 2L, n_p10_genes >= 3 ~ 3L, TRUE ~ 4L)) %>%
  filter(n_sig_genes >= 1 | n_p10_genes >= 2 | forced) %>%                                # <-- display rule
  filter(!str_detect(Species, regex("Bacteroides[_ ]fragilis", ignore_case = TRUE))) %>%    # not shown in the submitted panel
  arrange(desc(forced), tier, desc(n_sig_genes), desc(n_p10_genes), best_p) %>%
  slice_head(n = 31)                                                                        # <-- row cap
cat("Species displayed:", nrow(sp_rank), "\n")
save_table(sp_rank, "Fig4C_displayed_species")

rho_mat <- with(cors, tapply(rho, list(Species, Gene), identity))[sp_rank$Species, host_genes]
p_mat   <- with(cors, tapply(P,   list(Species, Gene), identity))[sp_rank$Species, host_genes]
rho0 <- rho_mat; rho0[is.na(rho0)] <- 0
row_order <- rownames(rho0)[hclust(dist(rho0))$order]
col_order <- unlist(lapply(levels(host_block), function(b) {
  g <- host_genes[host_block == b]; if (length(g) <= 1) return(g)
  g[hclust(dist(t(rho0[, g])))$order]
}))
rho_mat <- rho_mat[row_order, col_order]; p_mat <- p_mat[row_order, col_order]
sp_rank <- sp_rank[match(row_order, sp_rank$Species), ]
host_stat <- host_stat[match(col_order, host_stat$Gene), ]
category_mat <- category_mat[, col_order, drop = FALSE]
## annotation rows in the order of the submitted panel (any other gene set is appended)
category_order <- c("PI3K–AKT signaling with NF-κB", "Lymphocyte-mediated immunity", "T-cell differentiation",
                    "Integrin cell-surface interactions", "Inflammatory leukocyte migration",
                    "Colorectal adenoma signature", "Epithelial–mesenchymal transition")
category_mat <- category_mat[c(intersect(category_order, rownames(category_mat)), setdiff(rownames(category_mat), category_order)), , drop = FALSE]
column_split <- host_block[col_order]


## ---- 6. Figure 4C ------------------------------------------------------------------------

tax <- taxonomy[match(row_order, taxonomy$Species), ]
tax$Phylum[is.na(tax$Phylum)] <- "Unclassified"; tax$Family[is.na(tax$Family)] <- "Unclassified"
singleton_ok <- c("Fusobacteriaceae", "Bifidobacteriaceae", "Oscillospiraceae", "Veillonellaceae")   # keep even if only one species
fam_n <- table(tax$Family)
tax$Family[tax$Family %in% names(fam_n)[fam_n == 1] & !tax$Family %in% singleton_ok] <- "Others"

phylum_cols <- c(Actinobacteria = "#C2D3EA", Actinomycetota = "#C2D3EA", Bacteroidota = "#F1CCAD", Bacteroidetes = "#F1CCAD",
                 Firmicutes = "#C8DFC0", Bacillota = "#C8DFC0", Fusobacteria = "#E8BBB8", Fusobacteriota = "#E8BBB8",
                 Proteobacteria = "#D6C7E7", Pseudomonadota = "#D6C7E7", Unclassified = "#EEEEEE")
family_cols <- c(Lachnospiraceae = "#7FA7C9", Bacteroidaceae = "#E5A36B", Streptococcaceae = "#8FBE88", Actinomycetaceae = "#B89ABD",
                 Fusobacteriaceae = "#D77C78", Oscillospiraceae = "#C0A27A", Bifidobacteriaceae = "#78B7B1", Veillonellaceae = "#A99BD2",
                 Prevotellaceae = "#C9BE72", Unclassified = "#F1F1F1", Others = "#E3E3E3")
phylum_levels <- intersect(names(phylum_cols), unique(tax$Phylum)); family_levels <- intersect(names(family_cols), unique(tax$Family))
stopifnot(all(tax$Phylum %in% phylum_levels), all(tax$Family %in% family_levels))

max_rho <- max(0.3, min(1, max(abs(rho_mat), na.rm = TRUE)))
lfc_lim <- max(0.5, quantile(abs(host_stat$LFC), 0.95))
col_rho <- colorRamp2(c(-max_rho, 0, max_rho), correlation_colors)
col_lfc <- colorRamp2(c(-lfc_lim, 0, lfc_lim), c(response_colors[["non-pCR"]], "#F7F7F7", response_colors[["pCR"]]))
col_p   <- colorRamp2(c(0, 1.5, 3), c("#F7FBFF", "#7AA6D1", "#2C5AA0"))
common <- list(cluster_rows = FALSE, cluster_columns = FALSE, column_split = column_split, column_title = NULL,
               column_gap = unit(1.5, "mm"), show_column_names = FALSE, rect_gp = gpar(col = "white", lwd = 0.3))

ht_cat <- do.call(Heatmap, c(list(category_mat, name = "Host category", col = c("Not included" = "#F2F2F2", "Included" = "#303030"),
                                  row_names_gp = gpar(fontsize = 8), height = unit(nrow(category_mat) * 4, "mm")), common))
ht_lfc <- do.call(Heatmap, c(list(matrix(host_stat$LFC, 1, dimnames = list("LFC", col_order)), name = "LFC", col = col_lfc,
                                  heatmap_legend_param = list(title = expression(Log[2] * " FC")),
                                  row_names_gp = gpar(fontsize = 8), height = unit(4.5, "mm")), common))
ht_p <- do.call(Heatmap, c(list(matrix(host_stat$neglog10P, 1, dimnames = list("p-value", col_order)), name = "-log10 p", col = col_p,
                                heatmap_legend_param = list(title = expression(-Log[10] * italic(p))),
                                row_names_gp = gpar(fontsize = 8, fontface = "italic"), height = unit(4.5, "mm"),
                                cell_fun = function(j, i, x, y, w, h, fill) {
                                  s <- host_stat$stars[j]; if (nzchar(s)) grid.text(s, x, y, gp = gpar(fontsize = 7, fontface = "bold",
                                                                                                          col = ifelse(host_stat$neglog10P[j] >= 1.5, "white", "black")))
                                }), common))
right <- rowAnnotation(Phylum = factor(tax$Phylum, levels = phylum_levels), Family = factor(tax$Family, levels = family_levels),
                       col = list(Phylum = phylum_cols[phylum_levels], Family = family_cols[family_levels]),
                       simple_anno_size = unit(4.5, "mm"), gap = unit(0.8, "mm"), annotation_name_rot = 90,
                       annotation_name_gp = gpar(fontsize = 8), gp = gpar(col = "white", lwd = 0.3))
ht_main <- do.call(Heatmap, c(list(rho_mat, name = "Cor.", col = col_rho, na_col = "#F2F2F2", right_annotation = right,
                                   row_labels = pretty_species(rownames(rho_mat)), row_names_gp = gpar(fontsize = 8, fontface = "italic"),
                                   column_names_gp = gpar(fontsize = 8), column_names_rot = 90,
                                   width = unit(ncol(rho_mat) * 3.7, "mm"), height = unit(nrow(rho_mat) * 4.6, "mm"),
                                   layer_fun = function(j, i, x, y, w, h, fill) {
                                     s <- p_stars(p_mat[cbind(i, j)]); k <- nzchar(s)
                                     if (any(k)) grid.text(s[k], x[k], y[k], gp = gpar(fontsize = 7, fontface = "bold",
                                                                                      col = ifelse(abs(rho_mat[cbind(i, j)][k]) >= 0.65, "white", "black")))
                                   }), common[setdiff(names(common), "show_column_names")], list(show_column_names = TRUE)))
lgd_stars <- Legend(title = expression(italic(p) * "-value"), labels = c("p < 0.05", "p < 0.01", "p < 0.001"),
                    graphics = lapply(c("*", "**", "***"), function(s) function(x, y, w, h) grid.text(s, x, y, gp = gpar(fontsize = 9, fontface = "bold"))))

dev <- open_panel("Fig4C_species_host_correlations", width = 7.4, height = 10.5)
draw_panel(dev, draw(ht_cat %v% ht_lfc %v% ht_p %v% ht_main, heatmap_legend_side = "right", annotation_legend_side = "right",
                     annotation_legend_list = list(lgd_stars), merge_legends = TRUE, padding = unit(c(1, 2, 1, 2), "mm")))

message("Done: 13_Fig4C_Species_host_correlations.R")
