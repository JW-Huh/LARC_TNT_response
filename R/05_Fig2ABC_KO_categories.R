## =============================================================================
##  05_Fig2ABC_KO_categories.R
##  ---------------------------------------------------------------------------
##  PANELS      Figure 2A  KO-level response associations, one strip per
##                         functional category (ordered by P within category)
##              Figure 2B  category-level "robust perturbation score" versus
##                         size-matched random KO sets
##              Figure 2C  heatmap of the KOs with P < 0.05
##
##  PURPOSE     Ask whether response-associated microbial genes concentrate in
##              particular functional categories, using a rule-based KEGG
##              Orthology (KO) -> 19-category assignment (Supplementary Table 3).
##
##  INPUT       input/microbiome/humann_ko_cpm.tsv               HUMAnN 3 KO
##                    abundance (copies per million, unstratified), 42 samples
##              input/microbiome/ko_functional_categories.csv    frozen KO ->
##                    category assignment (2,228 rows = 2,127 unique KOs; 101 KOs
##                    belong to two categories), upper (5) and main (19) category
##              input/metadata/stool_samples.csv, patients.csv
##
##  OUTPUT      output/figures/Fig2A_KO_category_strips.{svg,png}
##              output/figures/Fig2B_KO_category_enrichment.{svg,png}
##              output/figures/Fig2C_KO_heatmap.{svg,png}
##              output/tables/Fig2A_KO_wilcoxon.csv            per-assignment statistics
##              output/tables/Fig2B_KO_category_enrichment.csv
##              output/cache/ko_baseline.RData                  (used by later scripts)
##
##  METHODS
##   * KO filter: KOs detected (CPM > 0) in >= 50 % of the 26 baseline samples
##     (3,509 KOs); of these, the 2,127 KOs with a category assignment form the
##     analysis background.
##   * KO test: two-sided Wilcoxon rank-sum test pCR vs non-pCR on CPM;
##     log2FC = log2(mean_pCR / mean_nonpCR).  stats::wilcox.test with its
##     default `exact = NULL`: exact when there are no ties, otherwise the
##     normal approximation with continuity correction.  Supplementary Table 2E
##     was generated with the exact *conditional* test (ties handled by the
##     permutation distribution of the rank sum, as in coin::wilcox_test with
##     distribution = "exact"); the two agree exactly for KOs without ties and
##     differ slightly for the few tied KOs (e.g. K02406 flagellin: 0.0069 here
##     vs 0.0054 in the table).  Neither the 2C KO set nor its stars change.
##   * Category score (2B): for a category with n KOs, the fraction of its KOs
##     with P below each of eight thresholds (0.005 ... 0.5) is averaged
##     (= mean empirical CDF of the P values at those thresholds).  The same
##     score is computed for 10,000 random draws of n assignment rows from the
##     2,228-row background (seed 123).  Reported: Z = (observed - mean_null) /
##     sd_null ("robust perturbation score") and the empirical one-sided P
##     = (1 + #null >= observed) / (10,001).  Categories with < 15 KOs are not
##     tested; categories with < 20 KOs are flagged (triangles).  Because the
##     background is sampled by assignment row, the 101 doubly assigned KOs
##     carry double weight - this is how the submitted analysis was run.
##   * Heatmap (2C): KOs with P < 0.05 (duplicates collapsed), row z-scores of
##     CPM, rows clustered by Manhattan distance / Ward.D2, samples ordered
##     within each response group by Euclidean / complete clustering.
##     Significance stars in the P-value annotation: *** P < 0.005, ** P < 0.01
##     (the thresholds of the submitted panel and its legend).
##
##  R PACKAGES  ggplot2, ComplexHeatmap 2.18.0, circlize, dplyr
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
})
report_versions(c("ggplot2", "dplyr", "ComplexHeatmap", "circlize"))
set_heatmap_options()                                   # legend font sizes (R/_common.R)


## ---- 1. Baseline KO abundance and the frozen category assignment ---------------

baseline <- read_stool_samples(timepoint = "Before")                           # 26 patients
ko_all <- read_feature_table("humann_ko_cpm.tsv", "KO", samples = baseline$SampleID)   # KOs x 26 (CPM)

prevalence_cut <- 0.50                                                         # <-- KO detection filter
ko_cpm <- ko_all[rowMeans(ko_all > 0) >= prevalence_cut, ]
cat("KOs detected in >=", 100 * prevalence_cut, "% of baseline samples:", nrow(ko_cpm), "\n")

annotation <- read.csv(input_file("microbiome", "ko_functional_categories.csv"), stringsAsFactors = FALSE)
category_levels <- unique(annotation$assigned_category[order(annotation$assigned_category_order)])
upper_levels    <- unique(annotation$upper_category[order(annotation$upper_category_order)])
annotation$assigned_category <- factor(annotation$assigned_category, levels = category_levels)
annotation$upper_category    <- factor(annotation$upper_category, levels = upper_levels)
stopifnot(all(annotation$KO %in% rownames(ko_cpm)), !anyDuplicated(paste(annotation$KO, annotation$assigned_category)))
cat("Category assignments:", nrow(annotation), "rows,", length(unique(annotation$KO)), "unique KOs\n")

## Colours of the five upper categories (as submitted)
upper_colors <- c("Ecological interaction" = "#71BE81", "Biomolecule metabolism" = "#DB9261",
                  "Energy metabolism" = "#E6B750", "Information processing" = "#7BAADD",
                  "Cellular processes" = "#A991CF")


## ---- 2. Per-KO Wilcoxon tests ----------------------------------------------------

is_pcr <- baseline$Response == "pCR"
ko_stats <- data.frame(
  KO = rownames(ko_cpm),
  P = apply(ko_cpm, 1, function(x) wilcox_p(x, baseline$Response)),
  Mean_pCR = rowMeans(ko_cpm[, is_pcr]), Mean_nonpCR = rowMeans(ko_cpm[, !is_pcr]),
  row.names = NULL
)
ko_stats$Log2FC <- log2(ko_stats$Mean_pCR / ko_stats$Mean_nonpCR)
stopifnot(all(is.finite(ko_stats$P)))

## one row per (KO, category) assignment
ko_results <- annotation %>%
  select(KO, KO_name, upper_category, assigned_category) %>%
  left_join(ko_stats, by = "KO")
save_table(ko_results, "Fig2A_KO_wilcoxon")

save(ko_cpm, ko_stats, baseline, file = cache_file("ko_baseline.RData"))


## ---- 3. Figure 2A: P-value strips per category ------------------------------------
## Each category is one strip; KOs are ordered by increasing P from left to
## right, and the strip width is proportional to the number of KOs.  The
## label above each strip is "significant / total" at P < 0.05.  Category
## names are placed below (two staggered rows) or above the strips, as in the
## submitted panel; edit `label_slot` to move a label.

ko_results <- ko_results %>%
  arrange(assigned_category, P) %>%
  group_by(assigned_category) %>%
  mutate(Rank = row_number(), neg_log10_P = -log10(P)) %>%
  ungroup()

strip_gap <- 8                                   # <-- empty KO-widths between strips
strips <- ko_results %>%
  group_by(assigned_category, upper_category) %>%
  summarise(n = n(), n_sig = sum(P < 0.05), .groups = "drop") %>%
  arrange(assigned_category) %>%
  mutate(x_start = cumsum(c(0, head(n + strip_gap, -1))), x_end = x_start + n, x_mid = (x_start + x_end) / 2,
         count_label = paste0(n_sig, "/", n))
print(strips, n = Inf)
ko_results$x <- strips$x_start[match(ko_results$assigned_category, strips$assigned_category)] + ko_results$Rank

## where each category name goes: -1 / -2 = first / second row below, +1 = above
label_slot <- c(
  "Quorum sensing/biofilm/colonization" = -1, "Motility/adhesion" = -2, "Xenobiotic metabolism/AMR" = 1,
  "Virulence/Secretion systems" = -2, "Antimicrobial peptide" = 1, "Amino acid" = -1,
  "Carbohydrate/glycan" = -2, "Nucleotide" = -1, "Cofactor/vitamin" = -2, "Lipid" = -1, "Energy" = -1,
  "Fiber fermentation and SCFA" = -2, "Energy-coupled ion transport" = 1,
  "Transcription/Translation/RNA processing" = -1, "DNA replication/repair" = -2,
  "Signal transduction/regulation" = -1, "General transport" = -1,
  "Cell envelope/structure/stress response" = -2, "Ion and metal homeostasis" = -1
)
wrap_name <- function(x, width = 18) vapply(x, function(z) paste(strwrap(gsub("/", "/ ", z), width), collapse = "\n"), character(1))
y_max <- max(ko_results$neg_log10_P)
strips <- strips %>% mutate(
  slot = label_slot[as.character(assigned_category)],
  label = gsub("/ ", "/", wrap_name(as.character(assigned_category))),
  y_label = case_when(slot == -1 ~ -0.16 * y_max, slot == -2 ~ -0.42 * y_max, TRUE ~ 1.28 * y_max),
  y_tick_end = ifelse(slot > 0, 1.12 * y_max, ifelse(slot == -2, -0.30 * y_max, -0.04 * y_max))
)

p_2a <- ggplot(ko_results, aes(x, neg_log10_P, colour = upper_category)) +
  ## strip separators
  geom_segment(data = strips, aes(x = x_start - strip_gap / 2, xend = x_start - strip_gap / 2, y = 0, yend = 1.12 * y_max),
               inherit.aes = FALSE, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  annotate("text", x = 1, y = -log10(0.05), label = "italic(p) == 0.05", parse = TRUE, hjust = 0, vjust = -0.4,
           size = 2.6, colour = "grey40") +
  geom_point(size = 0.7) +
  ## "significant / total" above every strip
  geom_text(data = strips, aes(x = x_mid, y = 1.06 * y_max, label = count_label), inherit.aes = FALSE, size = 2.5) +
  ## category names with a short connector to the strip
  geom_segment(data = filter(strips, slot != -1), aes(x = x_mid, xend = x_mid, y = ifelse(slot > 0, 1.12 * y_max, 0), yend = y_tick_end),
               inherit.aes = FALSE, colour = "grey50", linewidth = 0.3) +
  geom_text(data = strips, aes(x = x_mid, y = y_label, label = label), inherit.aes = FALSE, size = 2.4, lineheight = 0.85) +
  scale_colour_manual(values = upper_colors, name = "Category") +
  scale_x_continuous(expand = expansion(add = c(strip_gap / 2, 6))) +      # panel edge = first separator
  scale_y_continuous(breaks = 0:3, expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.55 * y_max, 1.4 * y_max), clip = "off") +
  labs(x = NULL, y = expression(-Log[10](italic(p)))) +
  theme_tnt(base_size = 9, legend = "top") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line = element_blank(),
        legend.key.size = unit(0.4, "lines"), legend.margin = margin(0, 0, 0, 0))
save_panel(p_2a, "Fig2A_KO_category_strips", width = 12, height = 3.6)

## Note: the y axis shows the -log10(P) range of the points only; the label rows
## below/above the strips lie outside that range (coord_cartesian with clip off).


## ---- 4. Figure 2B: competitive category enrichment ---------------------------------

thresholds     <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50)   # <-- P-value grid of the score
n_permutations <- 10000L                                               # <-- random sets per category
min_ko_tested  <- 15L                                                  #     categories below are skipped
low_count_ko   <- 20L                                                  #     categories below are flagged

## Score = mean over thresholds of the fraction of KOs with P < threshold
cdf_score <- function(p_sorted, n) mean(findInterval(thresholds, p_sorted) / n)

## Background: P value of every assignment row, in the annotation's row order
## (the order matters for exact reproduction of the seeded permutations)
p_background <- ko_stats$P[match(annotation$KO, ko_stats$KO)]

set.seed(123)                                                          # <-- permutation seed
enrichment <- list()
for (category in sort(levels(annotation$assigned_category))) {        # alphabetical iteration, as submitted
  observed_p <- sort(p_background[annotation$assigned_category == category])
  n_ko <- length(observed_p)
  if (n_ko < min_ko_tested) next

  observed <- cdf_score(observed_p, n_ko)
  null_scores <- replicate(n_permutations, cdf_score(sort(sample(p_background, n_ko)), n_ko))

  enrichment[[category]] <- data.frame(
    Category = category, N_KO = n_ko, Score = observed,
    Null_mean = mean(null_scores), Null_sd = sd(null_scores),
    Z = (observed - mean(null_scores)) / sd(null_scores),
    P = (1 + sum(null_scores >= observed)) / (n_permutations + 1)
  )
}
enrichment <- do.call(rbind, enrichment)
enrichment$Evidence <- ifelse(enrichment$N_KO < low_count_ko, "low-count category", "distributed signal")
enrichment <- enrichment[order(-enrichment$Z), ]
print(enrichment[, c("Category", "N_KO", "Z", "P", "Evidence")], row.names = FALSE, digits = 4)
save_table(enrichment, "Fig2B_KO_category_enrichment")

## Design (as submitted): lollipop of Z, circles (triangles for < 20 KOs) sized
## by KO count and filled by -log10(P) on a white-to-red gradient.
enrichment$Category <- factor(enrichment$Category, levels = rev(enrichment$Category))
p_2b <- ggplot(enrichment, aes(Z, Category)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = Z, yend = Category), colour = "grey75", linewidth = 0.4) +
  geom_point(aes(size = N_KO, fill = -log10(P), shape = Evidence), colour = "grey30", stroke = 0.4) +
  scale_shape_manual(values = c("distributed signal" = 21, "low-count category" = 24), name = "Evidence type") +
  scale_fill_gradient(low = "#F7E5E7", high = "#A74253", name = expression(-Log[10](italic(p)))) +
  scale_size_continuous(range = c(2, 7), breaks = c(100, 200, 300), name = "Number\nof KOs") +
  scale_y_discrete(labels = function(z) gsub("/ ", "/", wrap_name(z, 26), fixed = TRUE)) +   # wrap long names (helper from 2A)
  guides(fill = guide_colourbar(order = 1, barheight = unit(2.2, "cm"), barwidth = unit(0.35, "cm")),
         size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey30")),
         shape = guide_legend(order = 3, override.aes = list(size = 2.5))) +
  labs(x = "Robust perturbation score", y = NULL) +
  theme_tnt(base_size = 9, legend = "right") +
  theme(legend.box = "vertical", axis.text.y = element_text(lineheight = 0.8, size = 7.5),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7.5))
save_panel(p_2b, "Fig2B_KO_category_enrichment", width = 4.0, height = 5.2)


## ---- 5. Figure 2C: heatmap of response-associated KOs -------------------------------

heat_p_cut <- 0.05                                                     # <-- KO inclusion threshold
star_cuts  <- c(0.005, 0.01); star_symbols <- c("***", "**")          # <-- star thresholds (*** P < 0.005, ** P < 0.01)
selected <- ko_results %>% filter(P < heat_p_cut) %>% distinct(KO, .keep_all = TRUE)
cat("KOs with P <", heat_p_cut, ":", nrow(selected), "\n")
cat("Stars:", paste(sprintf("%s %s", selected$KO, p_stars(selected$P, star_cuts, star_symbols))[selected$P < max(star_cuts)], collapse = ", "), "\n")

z <- t(scale(t(ko_cpm[selected$KO, ])))                                # row z-scores of CPM
rownames(z) <- paste0(selected$KO, ": ", sub(" \\[EC:.*\\]$", "", selected$KO_name))   # "K03736: ethanolamine ammonia-lyase ..." 
stopifnot(all(is.finite(z)))

col_order <- order_samples_within_response(z, baseline)
meta <- baseline[match(col_order, baseline$SampleID), ]
z <- z[, col_order]

lfc <- selected$Log2FC
lfc[!is.finite(lfc)] <- sign(lfc[!is.finite(lfc)]) * max(abs(lfc[is.finite(lfc)]))

right <- rowAnnotation(
  `p-value` = anno_simple(-log10(selected$P), col = colorRamp2(c(1, 2.5), c("white", "#2166AC")),
                          pch = p_stars(selected$P, star_cuts, star_symbols),
                          pt_gp = gpar(col = "white", fontsize = 8), gp = gpar(col = "grey40", lwd = 0.5),
                          width = unit(3.5, "mm")),
  LFC = anno_simple(lfc, col = colorRamp2(c(-2, 0, 4), c("#D1606E", "white", "#57AD93")),
                    gp = gpar(col = "grey40", lwd = 0.5), width = unit(3.5, "mm")),
  Category = as.character(selected$upper_category),
  col = list(Category = upper_colors),
  gp = gpar(col = "grey40", lwd = 0.5), simple_anno_size = unit(3.5, "mm"),
  annotation_name_gp = gpar(fontsize = 8), annotation_name_rot = 90,
  annotation_legend_param = list(Category = list(ncol = 2))
)
ht <- Heatmap(
  z, name = "z-score", col = colorRamp2(c(-2, 0, 2), zscore_colors),
  cluster_columns = FALSE, column_split = meta$Response, column_title = NULL, column_gap = unit(1, "mm"),
  clustering_distance_rows = "manhattan", clustering_method_rows = "ward.D2",
  row_dend_side = "left", row_dend_width = unit(1.2, "cm"),
  show_column_names = FALSE, row_names_side = "right", row_names_gp = gpar(fontsize = 7),
  row_names_max_width = unit(12, "cm"), rect_gp = gpar(col = "grey40", lwd = 0.5),
  top_annotation = response_top_annotation(meta), right_annotation = right,
  heatmap_legend_param = list(direction = "horizontal", at = -2:2, title_position = "topcenter")
)
lgd_p   <- Legend(title = expression(-Log[10](italic(p))), col_fun = colorRamp2(c(1, 2.5), c("white", "#2166AC")),
                  at = c(1, 1.5, 2, 2.5), direction = "horizontal", title_position = "topcenter")
lgd_lfc <- Legend(title = expression(Log[2]*FC), col_fun = colorRamp2(c(-2, 0, 4), c("#D1606E", "white", "#57AD93")),
                  at = c(-2, 0, 2, 4), direction = "horizontal", title_position = "topcenter")

dev <- open_panel("Fig2C_KO_heatmap", width = 10, height = 8.5)
draw_panel(dev, draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                     merge_legends = TRUE, annotation_legend_list = list(lgd_p, lgd_lfc)))

message("Done: 05_Fig2ABC_KO_categories.R")
