## =============================================================================
##  11_Fig4B_Host_gene_set_enrichment.R
##  ---------------------------------------------------------------------------
##  PANEL       Figure 4B  gene-set enrichment (GSEA) of tumour expression,
##                         pCR versus non-pCR: 22 selected gene sets
##
##  PURPOSE     Summarise which host transcriptional programmes are shifted
##              between response groups, using pre-ranked GSEA on the DESeq2
##              Wald statistic and a curated catalogue of 2,736 MSigDB gene sets
##              grouped into host "axes" (immunity, barrier, ECM, ...).
##
##  TWO MODES   refit_gsea <- FALSE (default)  use the frozen GSEA result table
##                  that produced the submitted panel (fgsea is stochastic and
##                  its random state was not saved; the frozen table is the
##                  record).  The 22 displayed statistics are checked against
##                  the submitted values.
##              refit_gsea <- TRUE   rerun fgsea from output/cache/host_deseq2.RData
##                  (script 10; needs the restricted count data).  P values and
##                  FDR will differ slightly from the submission (a new run).
##
##  INPUT       input/host_rnaseq/gsea_targeted_results_frozen.csv   2,803 rows:
##                  gene set x host axis, fgsea statistics (pval, padj, ES, NES,
##                  size, leading edge)
##              input/host_rnaseq/msigdb_2026.1.Hs_gene_sets.rds   the gene sets
##                  (for a refit)
##              input/host_rnaseq/fig4C_gsea_annotation.csv        gene-set ->
##                  category annotation used by Figure 4C (passed on)
##
##  OUTPUT      output/figures/Fig4B_host_gene_sets.{svg,png}
##              output/tables/Fig4B_selected_gene_sets.csv
##              output/cache/host_enrichment.RData   (fgsea table, 4C annotation)
##
##  METHODS     fgsea 1.28.0 (multilevel algorithm, eps = 0, minSize 5, maxSize
##              500, seed 123) on genes ranked by the DESeq2 Wald statistic
##              (positive = higher in pCR).  NES = enrichment score normalised
##              to the mean of the permutation null; FDR = Benjamini-Hochberg
##              over the whole 2,803-row testing family (not over the 22 shown).
##              The 22 gene sets are a curated selection of significant sets
##              (FDR < 0.1, |NES| >= 1) representing distinct programmes.
##
##  R PACKAGES  fgsea (refit only), dplyr, ggplot2, patchwork
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
})
report_versions(c("fgsea", "ggplot2", "patchwork"))

refit_gsea <- FALSE        # <-- TRUE = rerun fgsea from the DESeq2 ranks (needs script 10)


## ---- 1. GSEA result table --------------------------------------------------------

if (!refit_gsea) {
  fgsea_res <- read.csv(input_file("host_rnaseq", "gsea_targeted_results_frozen.csv"), stringsAsFactors = FALSE)
  message("Using the frozen GSEA table (", nrow(fgsea_res), " gene set x axis rows).")
} else {
  suppressPackageStartupMessages(library(fgsea))
  load(cache_file("host_deseq2.RData"))                                   # deg (script 10)
  gs <- readRDS(input_file("host_rnaseq", "msigdb_2026.1.Hs_gene_sets.rds"))
  ranks <- deg %>% filter(is.finite(stat)) %>% arrange(desc(abs(stat))) %>% distinct(Gene, .keep_all = TRUE)
  ranks <- setNames(ranks$stat, ranks$Gene)
  set.seed(123)
  fgsea_res <- fgsea(pathways = gs$gene_sets, stats = ranks, minSize = 5, maxSize = 500, eps = 0, nproc = 1) %>%
    as.data.frame() %>%
    mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)), pathway_id = pathway) %>%
    left_join(gs$gene_set_table, by = "pathway_id") %>%          # adds Host_axis, gs_name, description, sizes
    arrange(padj, pval, desc(abs(NES)))
  message("fgsea rerun: ", nrow(fgsea_res), " rows (a NEW run; statistics differ from the submission).")
}
stopifnot(all(c("gs_name", "NES", "pval", "padj") %in% names(fgsea_res)))


## ---- 2. The 22 gene sets of Figure 4B ---------------------------------------------
## `gs_name` is the exact MSigDB identifier that was tested; `display_label` is
## the wording of the submitted panel.  NOTE: the first row's label in the
## submission ("PI3K-AKT signaling with NF-kB") refers to the KEGG MEDICUS set
## "NNK/NNN -> PI3K signaling (N01339)"; change `display_label` if you prefer
## the literal set name.  `expected_*` are the submitted statistics used as a
## consistency check in replay mode.

catalog <- tibble::tribble(
  ~gs_name, ~display_label, ~main_category, ~expected_NES, ~expected_P, ~expected_FDR,
  "KEGG_MEDICUS_ENV_FACTOR_NNK_NNN_TO_PI3K_SIGNALING_PATHWAY_N01339", "PI3K–AKT signaling with NF-κB", "Growth / vascular signaling", 1.98714290212083, 0.000514610429646659, 0.0198306216269473,
  "GOMF_NAD_P_H_OXIDASE_H2O2_FORMING_ACTIVITY", "NAD(P)H oxidase activity", "Redox / metabolic signaling", 1.93090891572615, 0.00123415887408428, 0.0312653581434684,
  "REACTOME_O_LINKED_GLYCOSYLATION_OF_MUCINS", "Mucin glycosylation", "Barrier / mucin", 1.84090549547323, 0.000492340503964239, 0.0195663948271799,
  "HALLMARK_G2M_CHECKPOINT", "G2/M checkpoint", "Cell cycle / tumor signatures", 1.75348890339648, 2.11184286446835e-05, 0.00240750086549392,
  "GOBP_ALPHA_BETA_T_CELL_ACTIVATION", "αβ T-cell activation", "Adaptive immunity", 1.70926732667661, 6.91690928768895e-05, 0.00540704680317628,
  "GOBP_T_CELL_DIFFERENTIATION", "T-cell differentiation", "Adaptive immunity", 1.69554999047592, 8.03479359922496e-06, 0.00118237892010027,
  "WP_ARYL_HYDROCARBON_RECEPTOR_PATHWAY_WP2586", "AhR pathway", "Redox / metabolic signaling", 1.6535695451036, 0.00443756195227686, 0.0663120074781133,
  "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY", "Lymphocyte-mediated immunity", "Adaptive immunity", 1.6320836266624, 7.55800222640287e-06, 0.00118237892010027,
  "KEGG_TIGHT_JUNCTION", "Tight junction", "Barrier / mucin", 1.6161033146903, 0.00199674656796764, 0.0401698427202902,
  "GOBP_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE", "Innate immune activation", "Innate immunity", 1.55771938640985, 0.00026624860903863, 0.0134899295246239,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "IFN-γ response", "Innate immunity", 1.48816851075625, 0.00230332713027794, 0.044555129162717,
  "WP_VEGFAVEGFR2_SIGNALING", "VEGFA–VEGFR2 signaling", "Growth / vascular signaling", 1.41782696906224, 0.000849972390039476, 0.0255552138367913,
  "HALLMARK_MYC_TARGETS_V1", "MYC targets", "Cell cycle / tumor signatures", -1.50817303765219, 0.00169653563938335, 0.0357055500719449,
  "GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE", "Antimicrobial humoral response", "Innate immunity", -1.59521810369097, 0.00104489665512429, 0.0272270214135244,
  "NABA_CORE_MATRISOME", "Core matrisome", "ECM / EMT", -1.6196825579451, 0.000124155508718382, 0.00754865493007762,
  "SANSOM_WNT_PATHWAY_REQUIRE_MYC", "Wnt target signature", "Cell cycle / tumor signatures", -1.63054356140198, 0.00452746060365204, 0.0663120074781133,
  "WP_COMPLEMENT_SYSTEM", "Complement system", "Innate immunity", -1.67589050868522, 0.00106823008834064, 0.02757242945,
  "GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT", "Extracellular matrix structural support", "ECM / EMT", -1.74054990565147, 6.83669332630373e-05, 0.00540704680317628,
  "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS", "Integrin cell-surface interactions", "Cell adhesion / trafficking", -1.88635905362245, 9.40636638895285e-05, 0.00623434616425524,
  "GOBP_LEUKOCYTE_MIGRATION_INVOLVED_IN_INFLAMMATORY_RESPONSE", "Inflammatory leukocyte migration", "Cell adhesion / trafficking", -2.06427371888218, 0.000607213786299238, 0.0209620217514878,
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "Epithelial–mesenchymal transition", "ECM / EMT", -2.17974827269977, 2.52379206134686e-10, 1.381019015969e-07,
  "SABATES_COLORECTAL_ADENOMA_UP", "Colorectal adenoma-associated signature", "Cell cycle / tumor signatures", -2.54267740376741, 1.99576791126236e-13, 1.36510525130346e-10
)
stopifnot(nrow(catalog) == 22, !anyDuplicated(catalog$gs_name))

## one row per gene set (a set may be listed under several host axes with identical statistics)
selected <- fgsea_res %>%
  filter(gs_name %in% catalog$gs_name) %>%
  group_by(gs_name) %>%
  summarise(NES = NES[1], P = pval[1], FDR = padj[1], size = size[1],
            n_axes = n(), stat_conflict = max(abs(NES - NES[1]), abs(pval - pval[1])) > 1e-12, .groups = "drop")
stopifnot(setequal(selected$gs_name, catalog$gs_name), !any(selected$stat_conflict))
selected <- catalog %>% left_join(selected, by = "gs_name")

if (!refit_gsea) {
  dev <- with(selected, pmax(abs(NES - expected_NES), abs(P - expected_P), abs(FDR - expected_FDR)))
  stopifnot(all(dev < 1e-8))                 # frozen table == submitted statistics
  message("All 22 displayed statistics match the submitted panel.")
} else {
  message("Refit: displayed statistics are new; expected_* columns show the submitted values for comparison.")
}
stopifnot(all(selected$FDR < 0.10), all(abs(selected$NES) >= 1))
save_table(selected, "Fig4B_selected_gene_sets")


## ---- 3. Figure 4B -----------------------------------------------------------------
## Design (as submitted): gene-set names | category tile | FDR tile with stars |
## lollipop of NES coloured by direction; legends on the right.

category_order <- c("Adaptive immunity", "Innate immunity", "Barrier / mucin", "Cell adhesion / trafficking",
                    "Cell cycle / tumor signatures", "ECM / EMT", "Growth / vascular signaling", "Redox / metabolic signaling")
category_colors <- c("Adaptive immunity" = "#7A9CC6", "Innate immunity" = "#A8C5DF", "Barrier / mucin" = "#B6D4CA",
                     "Cell adhesion / trafficking" = "#70B7AA", "Cell cycle / tumor signatures" = "#DA8886",
                     "ECM / EMT" = "#E5ABA4", "Growth / vascular signaling" = "#C86F77", "Redox / metabolic signaling" = "#B5C98F")

d <- selected %>%
  arrange(NES) %>%                                                          # highest NES at the top
  mutate(y = row_number(), main_category = factor(main_category, levels = category_order),
         Direction = factor(ifelse(NES > 0, "pCR", "non-pCR"), levels = response_levels),
         logFDR = pmin(-log10(FDR), 10), stars = p_stars(FDR),
         star_colour = ifelse(logFDR >= 5, "white", "grey15"))

y_scale <- function() scale_y_continuous(limits = c(0.5, nrow(d) + 0.5), expand = c(0, 0))
tile_theme <- theme_void(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5), plot.margin = margin(3, 1, 3, 1))

p_labels <- ggplot(d, aes(1, y, label = display_label)) +
  geom_text(hjust = 1, size = 2.9) +
  y_scale() + scale_x_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  theme_void(base_size = 8) + theme(plot.margin = margin(3, 0, 3, 0))

p_category <- ggplot(d, aes(1, y, fill = main_category)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  y_scale() + scale_x_continuous(breaks = 1, labels = "Category", expand = c(0, 0)) +
  scale_fill_manual(values = category_colors, limits = category_order, drop = FALSE, name = "Category") +
  guides(fill = guide_legend(order = 1)) + tile_theme

p_fdr <- ggplot(d, aes(1, y, fill = logFDR)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = stars, colour = star_colour), size = 2) +
  scale_colour_identity() +
  y_scale() + scale_x_continuous(breaks = 1, labels = "FDR", expand = c(0, 0)) +
  scale_fill_gradient(low = "#F4F4F8", high = "#08519C", limits = c(0, 10), breaks = c(0, 2.5, 5, 7.5, 10),
                      name = expression(-Log[10] ~ FDR)) +
  guides(fill = guide_colourbar(order = 2, barwidth = unit(3.5, "mm"), barheight = unit(26, "mm"))) + tile_theme

p_nes <- ggplot(d, aes(y = y, colour = Direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_segment(aes(x = 0, xend = NES, yend = y), linewidth = 0.5, alpha = 0.7) +
  geom_point(aes(x = NES), size = 3) +
  y_scale() + scale_x_continuous(limits = c(-2.95, 2.35), breaks = -2:2, expand = c(0, 0)) +
  scale_colour_manual(values = response_colors, breaks = response_levels, name = "Enriched in") +
  guides(colour = guide_legend(order = 3)) +
  labs(x = "NES", y = NULL) +
  theme_tnt(base_size = 8) +
  theme(axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = "grey25", fill = NA, linewidth = 0.5),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3), plot.margin = margin(3, 2, 3, 1))

p_4b <- wrap_plots(p_labels, p_category, p_fdr, p_nes, widths = c(5.5, 0.3, 0.3, 2.9), guides = "collect") &
  theme(legend.position = "right", legend.title = element_text(size = 9), legend.text = element_text(size = 8),
        legend.key.height = unit(3.7, "mm"), legend.key.width = unit(3.7, "mm"))
save_panel(p_4b, "Fig4B_host_gene_sets", width = 7.2, height = 5.6)


## ---- 4. Pass the gene-set -> category annotation on to Figure 4C ------------------

gsea_annotation <- read.csv(input_file("host_rnaseq", "fig4C_gsea_annotation.csv"), stringsAsFactors = FALSE)
hit <- match(gsea_annotation$gs_name, catalog$gs_name)
gsea_annotation$display_label[!is.na(hit)] <- catalog$display_label[hit[!is.na(hit)]]
save(fgsea_res, gsea_annotation, file = cache_file("host_enrichment.RData"))

message("Done: 11_Fig4B_Host_gene_set_enrichment.R")
