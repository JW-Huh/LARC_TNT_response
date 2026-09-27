## =============================================================================
##  15_Fig4D_Metabolite_host_network.R
##  ---------------------------------------------------------------------------
##  PANEL       Figure 4D  bipartite "rail" network: fecal metabolites (centre)
##                         linked to pCR-enriched (left) and non-pCR-enriched
##                         (right) tumour genes by Spearman correlation
##
##  PURPOSE     Show which response-associated host genes co-vary with fecal
##              metabolite levels across the 17 patients that have both data
##              types, and whether the association survives adjustment for
##              response status itself (TRG).
##
##  INPUT       output/cache/host_metabolite_inputs.RData   (script 14)
##
##  OUTPUT      output/figures/Fig4D_metabolite_host_network.{svg,png}
##              output/tables/Fig4D_metabolite_host_correlations.csv   all pairs tested
##              output/tables/Fig4D_displayed_edges.csv, Fig4D_displayed_nodes.csv
##
##  METHODS
##   * Features: metabolites passing the QC of script 14, assay floor (minimum)
##     subtracted; host genes with DESeq2 P < 0.1 or Wilcoxon P < 0.1.
##   * Correlation: Spearman rho (Pearson correlation of ranks), t-distribution
##     P, n >= 8; BH FDR over all pairs.
##   * TRG-adjusted correlation: ranks of both variables are regressed on the
##     response group (model.matrix ~ Group) and the residuals correlated
##     (df = n - 2 - 1).  "Support" = same sign as the unadjusted rho, |rho| >=
##     0.25, P < 0.05, and the same sign in >= 80 % of leave-one-out refits.
##   * Display selection (a curated, rule-based reduction of ~1,000 nominal
##     edges): edges with |rho| >= 0.3 and P < 0.1 to genes with a functional
##     annotation (regex rules below); metabolites ranked by differential tier,
##     at most 2 per chemical class and 8 in total (butyrate always kept,
##     lithocholic acid as literature priority, kynurenic acid excluded);
##     up to 24 pCR- and 20 non-pCR-enriched genes; edges added greedily by P
##     with at most 8 edges per metabolite and 2 per gene; metabolites with
##     fewer than 3 displayed genes dropped (except butyrate).
##
##  R PACKAGES  dplyr, tidyr, stringr, ggplot2, ggnewscale
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggnewscale)
})
report_versions(c("ggplot2", "ggnewscale", "dplyr"))

load(cache_file("host_metabolite_inputs.RData"))          # meta, met_raw, host_expr, host_deg ...
stopifnot(identical(rownames(met_raw), meta$RNA_sample_id), identical(rownames(host_expr), meta$RNA_sample_id))

## ---- 0. Rules and thresholds (edit here) -----------------------------------------------

gene_p_cut        <- 0.10      # host gene enters if DESeq2 P or Wilcoxon P below this
edge_p_cut        <- 0.10      # nominal Spearman P for an edge
edge_rho_cut      <- 0.30      # |rho| for an edge
adj_p_cut         <- 0.05      # TRG-adjusted support: P
adj_rho_cut       <- 0.25      #                       |rho|
loo_sign_cut      <- 0.80      #                       leave-one-out sign agreement
min_n             <- 8L
max_metabolites   <- 8L;  max_per_class <- 2L
max_genes         <- 44L; max_edges_per_metabolite <- 8L; max_edges_per_gene <- 2L
priority_metabolites   <- "butyrate"                  # always displayed
literature_metabolites <- "lithocholic acid"          # displayed if it has any qualifying edge
excluded_metabolites   <- "kynurenic acid"

priority_genes <- c("CCR1", "TJP1", "S100A9", "CXCL1", "MARCO", "SIPA1", "BATF3", "LGR5", "GADD45B", "BOK",
                    "RUNX3", "NSD3", "GPR137B", "SART3", paste0("IGFBP", 1:7))
excluded_genes <- c("SENP3-EIF4A1", "DGCR11", "SP2-AS1", "C9ORF16", "UCA1", "MALAT1", "NEAT1", "XIST")
noncoding_regex <- regex(paste0("-AS[0-9]*$|^LINC[0-9]+$|^LOC[0-9]+$|^C[0-9XY]+ORF[0-9]+$|^AC[0-9]+|^AL[0-9]+|^AP[0-9]+|",
                                "^RP[0-9]+-|^MIR[0-9]+|^SNOR[AD][0-9]+|^RNU[0-9]+|^MT-T|^MT-R|^HCG[0-9]+|",
                                "^RPL[0-9A-Z]*P[0-9]+$|^RPS[0-9A-Z]*P[0-9]+$|^KRT[0-9A-Z]*P[0-9]+$|^EEF1A1P[0-9]+$|",
                                "^GAPDHP[0-9]+$|^HNRNP[A-Z0-9]*P[0-9]+$"), ignore_case = TRUE)

category_colors <- c("Tumor remodeling/CRC" = "#DF8A82", "Treatment response" = "#F4A261", "Metabolite signaling" = "#A8C98B",
                     "Barrier/mucus/AMP" = "#D9B44A", "Innate immunity" = "#9DC7DD", "Adaptive immunity" = "#B7A3CE",
                     "Gut homing" = "#77C1B5")

## functional category of a gene: manual table first, then regex rules (first match wins)
manual_category <- tibble::tribble(
  ~Gene, ~Category, ~priority,
  "CCR1", "Gut homing", 4, "TJP1", "Barrier/mucus/AMP", 4, "S100A9", "Innate immunity", 4, "CXCL1", "Innate immunity", 4,
  "MARCO", "Innate immunity", 4, "SIPA1", "Adaptive immunity", 4, "BATF3", "Adaptive immunity", 4, "LGR5", "Tumor remodeling/CRC", 4,
  "GADD45B", "Treatment response", 4, "BOK", "Treatment response", 4, "RUNX3", "Adaptive immunity", 4, "NSD3", "Tumor remodeling/CRC", 3,
  "GPR137B", "Metabolite signaling", 3, "SART3", "Treatment response", 3
) %>% bind_rows(tibble(Gene = paste0("IGFBP", 1:7), Category = "Treatment response", priority = 4))
category_rules <- tibble::tribble(
  ~regex, ~Category, ~priority,
  "^(LGR5|AXIN2|NOTUM|NKD1|MYC|MKI67|TOP2A|CDK1|CCNB1|CCNB2|BUB1|BUB1B|AURKA|AURKB|CDC20|PLK1)$", "Tumor remodeling/CRC", 3,
  "^(SNAI1|SNAI2|ZEB1|ZEB2|TWIST1|VIM|FN1|MMP1|MMP2|MMP3|MMP7|MMP9|MMP14|COL1A1|COL1A2|COL3A1|COL22A1|COMP|PDLIM4|CRK|MSX1|NSD3)$", "Tumor remodeling/CRC", 3,
  "^(BAX|BAK1|BCL2|BCL2L1|BOK|CASP3|CASP7|CASP8|CASP9|FAS|FASLG|TNFAIP1|GADD45A|GADD45B|GADD45G|IGFBP[1-7])$", "Treatment response", 3,
  "^(SOD2|HMOX1|NQO1|TXNRD1|GPX2|NOX1|DUOX2|ATM|ATR|BRCA1|BRCA2|RAD51|PARP1)$", "Treatment response", 3,
  "^(AHR|ARNT|CYP1A1|CYP1B1|NR1H4|GPBAR1|NR0B2|FGF19|SLC10A2|CYP7A1|S1PR2|S1PR3|HCAR2|HCAR3|NIACR1|FFAR2|FFAR3|GPR41|GPR43|GPR137B)$", "Metabolite signaling", 3,
  "^(SLC2A1|HK2|LDHA|PDK1|CA9|VEGFA|ENO1|ALDOA|PKM|SLC16A1|SLC16A3|SLC16A4)$", "Metabolite signaling", 2,
  "^(TJP1|TJP2|TJP3|OCLN|EPCAM|SLC9A3|MUC1|MUC2|MUC5AC|MUC5B|REG1A|REG1B|REG3A|REG3G|DEFA[1-9]|DEFB[0-9]+|CLDN[1-9][0-9]*)$", "Barrier/mucus/AMP", 3,
  "^(S100A8|S100A9|MARCO|TLR[1-9][0-9]*|NOD1|NOD2|MYD88|IRAK[1-4]|NLRP3|IL1B|TNF|CXCL1|CXCL2|CXCL3|CXCL5|CXCL8|CCL2|CCL3|CCL4|CCL5|C3|C5AR1|FCGR[1-3][A-Z]*)$", "Innate immunity", 3,
  "^(BATF3|SIPA1|RUNX3|CD3D|CD3E|CD3G|CD4|CD8A|CD8B|TRAC|TRBC1|TRBC2|KLRK1|NKG7|GNLY|GZMB|PRF1|TBX21|GATA3|FOXP3|IL7R)$", "Adaptive immunity", 3,
  "^(CCR1|CCR2|CCR5|CCR6|CCR7|CCR9|CXCR3|CXCR4|CX3CR1|ITGA4|ITGB7|ICAM1|VCAM1|SELE|SELL|SELPLG|CCL25|CX3CL1)$", "Gut homing", 3
)
gene_category <- function(genes) {
  cat <- manual_category$Category[match(genes, manual_category$Gene)]
  pri <- manual_category$priority[match(genes, manual_category$Gene)]
  for (i in seq_len(nrow(category_rules))) {
    hit <- is.na(cat) & str_detect(genes, regex(category_rules$regex[i], ignore_case = TRUE))
    cat[hit] <- category_rules$Category[i]; pri[hit] <- category_rules$priority[i]
  }
  data.frame(Gene = genes, Category = cat, priority = ifelse(is.na(pri), 0, pri))
}

metabolite_class <- function(key) case_when(
  key %in% c("acetate", "propionate", "butyrate", "isobutyrate", "isovalerate", "valerate", "gamma aminobutyric acid") ~ "SCFA / related",
  str_detect(key, "cholic|deoxycholic|ursodeoxycholic|lithocholic|tauro|glyco") ~ "Bile acid",
  key %in% c("indolepropionic acid", "indole lactic acid", "indole acetic acid", "indole", "tryptamine", "tryptophan",
             "kynurenic acid", "xanthurenic acid", "nicotinic acid") ~ "Tryptophan / indole",
  TRUE ~ "Other metabolite")


## ---- 1. Feature tables ------------------------------------------------------------------

## host genes: differential-expression evidence in this subset
genes <- host_deg %>%
  mutate(min_p = pmin(DESeq2_P, Wilcoxon_P, na.rm = TRUE)) %>%
  filter(is.finite(LFC), (is.finite(DESeq2_P) & DESeq2_P < gene_p_cut) | (is.finite(Wilcoxon_P) & Wilcoxon_P < gene_p_cut),
         Gene %in% colnames(host_expr)) %>%
  arrange(min_p, desc(abs(LFC))) %>% distinct(Gene, .keep_all = TRUE)
genes <- left_join(genes, gene_category(genes$Gene), by = "Gene")        # functional category (may be NA)
expr <- host_expr[, genes$Gene, drop = FALSE]

## metabolites: readable keys, assay floor subtracted, non-constant
met <- met_raw
colnames(met) <- str_to_lower(str_squish(gsub("_", " ", colnames(met))))
floor <- apply(met, 2, min)
met <- sweep(met, 2, ifelse(floor > 0, floor, 0))                     # subtract the assay floor
keep <- apply(met, 2, function(x) mean(!is.finite(x)) <= 0.2 & sum(is.finite(x)) >= min_n & length(unique(x)) > 1)
met <- met[, keep, drop = FALSE]
cat("Features:", ncol(met), "metabolites x", ncol(expr), "host genes,", nrow(meta), "patients\n")

is_pcr <- meta$Response == "pCR"
met_nodes <- data.frame(Metabolite = colnames(met),
                        log2FC = apply(met, 2, function(x) log2((median(x[is_pcr]) + 1e-6) / (median(x[!is_pcr]) + 1e-6))),
                        P = apply(met, 2, function(x) wilcox_p(x, meta$Response, exact = FALSE)), row.names = NULL) %>%
  mutate(FDR = p.adjust(P, "BH"), Class = metabolite_class(Metabolite), Label = str_to_sentence(Metabolite))


## ---- 2. Spearman correlations and TRG-adjusted correlations ------------------------------

## Pearson correlation of two matrices with t-test P values
cor_matrices <- function(x, y, df_loss = 0) {
  rho <- suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
  n <- crossprod(1 * is.finite(x), 1 * is.finite(y)); df <- n - 2 - df_loss
  p <- 2 * pt(-abs(rho * sqrt(pmax(df, 0) / pmax(1 - rho^2, .Machine$double.eps))), df = pmax(df, 1))
  bad <- n < min_n | !is.finite(rho) | df <= 0; rho[bad] <- NA; p[bad] <- NA
  list(rho = rho, p = p, n = n)
}
rank_cols <- function(m) apply(m, 2, rank, ties.method = "average")

raw <- cor_matrices(rank_cols(met), rank_cols(expr))                    # Spearman = Pearson on ranks

design <- model.matrix(~ Response, data = meta)                          # TRG adjustment
residualise <- function(m) apply(m, 2, function(x) lm.fit(design, rank(x, ties.method = "average"))$residuals)
met_res <- residualise(met); expr_res <- residualise(expr)
adj <- cor_matrices(met_res, expr_res, df_loss = qr(design)$rank - 1)

pairs <- expand_grid(Metabolite = colnames(met), Gene = colnames(expr)) %>%
  mutate(n = raw$n[cbind(Metabolite, Gene)], rho = raw$rho[cbind(Metabolite, Gene)], P = raw$p[cbind(Metabolite, Gene)],
         FDR = p.adjust(P, "BH"),
         adj_rho = adj$rho[cbind(Metabolite, Gene)], adj_P = adj$p[cbind(Metabolite, Gene)]) %>%
  mutate(adj_FDR = p.adjust(adj_P, "BH")) %>%
  filter(is.finite(rho))
save_table(pairs, "Fig4D_metabolite_host_correlations")
cat("Pairs tested:", nrow(pairs), "; nominal P < 0.1 & |rho| >= 0.3:", sum(pairs$P < edge_p_cut & abs(pairs$rho) >= edge_rho_cut), "\n")


## ---- 3. Display selection -------------------------------------------------------------

edges <- pairs %>%
  filter(abs(rho) >= edge_rho_cut, P < edge_p_cut) %>%
  left_join(genes %>% select(Gene, LFC, DESeq2_P, DESeq2_FDR, min_p, Category, priority), by = "Gene") %>%
  left_join(met_nodes, by = "Metabolite", suffix = c("", "_met")) %>%
  mutate(high_confidence = P < 0.01 | FDR < 0.10,
         excluded_gene = Gene %in% excluded_genes | str_detect(Gene, noncoding_regex),
         priority_gene = Gene %in% priority_genes,
         adj_support = is.finite(adj_rho) & sign(adj_rho) == sign(rho) & abs(adj_rho) >= adj_rho_cut & adj_P < adj_p_cut) %>%
  filter(!excluded_gene, !is.na(Category), priority > 0)             # annotated genes only

## 3a. metabolites
met_sel <- edges %>%
  group_by(Metabolite, Class, log2FC, P_met, FDR_met, Label) %>%
  summarise(n_genes = n_distinct(Gene), n_priority_genes = n_distinct(Gene[priority_gene]), sum_priority = sum(priority),
            n_high = sum(high_confidence), best_p = min(P), mean_abs_rho = mean(abs(rho)), .groups = "drop") %>%
  mutate(manual = Metabolite %in% priority_metabolites, literature = Metabolite %in% literature_metabolites,
         tier = case_when(FDR_met < 0.10 ~ 4L, P_met < 0.05 ~ 3L, P_met < 0.10 ~ 2L, abs(log2FC) >= 0.5 ~ 1L, TRUE ~ 0L)) %>%
  group_by(Class) %>%
  arrange(desc(tier), desc(abs(log2FC)), desc(n_priority_genes), desc(sum_priority), desc(n_high), best_p, .by_group = TRUE) %>%
  mutate(class_rank = row_number()) %>% ungroup() %>%
  filter(!Metabolite %in% excluded_metabolites, (n_genes >= 3 & class_rank <= max_per_class) | manual | literature) %>%
  arrange(desc(manual), desc(tier), desc(abs(log2FC)), desc(literature), desc(n_priority_genes), desc(sum_priority), desc(n_high), best_p, desc(mean_abs_rho)) %>%
  slice_head(n = max_metabolites) %>% mutate(selection_order = row_number())

## 3b. genes
pool <- edges %>% filter(Metabolite %in% met_sel$Metabolite)
gene_sel <- pool %>%
  group_by(Gene, LFC, DESeq2_FDR, Category, priority_gene) %>%
  summarise(priority = max(priority), n_met = n_distinct(Metabolite), n_high = sum(high_confidence), best_p = min(P),
            max_abs_rho = max(abs(rho)), .groups = "drop") %>%
  mutate(Side = ifelse(LFC > 0, "pCR", "non-pCR")) %>%
  arrange(best_p, desc(max_abs_rho), desc(n_high), desc(n_met), desc(priority), DESeq2_FDR, desc(priority_gene), Gene)
gene_sel <- bind_rows(filter(gene_sel, Side == "pCR") %>% slice_head(n = 24),
                      filter(gene_sel, Side == "non-pCR") %>% slice_head(n = 20), gene_sel) %>%
  distinct(Gene, .keep_all = TRUE) %>% slice_head(n = max_genes)

## 3c. edges: best edge of every gene, then greedy addition under the degree caps
pool <- pool %>% filter(Gene %in% gene_sel$Gene) %>% mutate(edge_id = paste(Metabolite, Gene, sep = "|"))
shown <- pool %>% group_by(Gene) %>% arrange(P, desc(abs(rho)), desc(high_confidence), desc(priority), desc(priority_gene), .by_group = TRUE) %>%
  slice_head(n = 1) %>% ungroup()
max_edges <- min(nrow(pool), nrow(met_sel) * max_edges_per_metabolite, nrow(gene_sel) * max_edges_per_gene)
while (nrow(shown) < max_edges) {
  candidates <- pool %>% filter(!edge_id %in% shown$edge_id) %>%
    left_join(count(shown, Metabolite, name = "deg_m"), by = "Metabolite") %>%
    left_join(count(shown, Gene, name = "deg_g"), by = "Gene") %>%
    filter(coalesce(deg_m, 0L) < max_edges_per_metabolite, coalesce(deg_g, 0L) < max_edges_per_gene) %>%
    arrange(P, desc(abs(rho)), desc(high_confidence), desc(priority), DESeq2_FDR, desc(priority_gene))
  if (nrow(candidates) == 0) break
  shown <- bind_rows(shown, candidates[1, names(shown)])
}
shown <- shown %>% group_by(Metabolite) %>% mutate(n_shown = n_distinct(Gene)) %>% ungroup() %>%
  filter(n_shown >= 3 | Metabolite %in% priority_metabolites)

met_sel <- met_sel %>% filter(Metabolite %in% shown$Metabolite) %>%
  arrange(desc(tier), desc(abs(log2FC)), selection_order) %>% mutate(order = row_number())
## submitted panel: glutamic acid was drawn just above indole acetic acid
if (all(c("glutamic acid", "indole acetic acid") %in% met_sel$Metabolite)) {
  met_sel$order[met_sel$Metabolite == "glutamic acid"] <- met_sel$order[met_sel$Metabolite == "indole acetic acid"] - 0.5
  met_sel <- arrange(met_sel, order)
}
gene_sel <- gene_sel %>% filter(Gene %in% shown$Gene) %>%
  left_join(shown %>% arrange(Gene, P, desc(abs(rho))) %>% group_by(Gene) %>% slice_head(n = 1) %>% ungroup() %>%
              transmute(Gene, anchor = Metabolite, anchor_order = match(Metabolite, met_sel$Metabolite), anchor_p = P, anchor_rho = abs(rho)),
            by = "Gene")

## 3d. leave-one-out sign agreement of the adjusted correlation for the displayed edges
shown$loo_sign <- mapply(function(m, g) {
  x <- met_res[, m]; y <- expr_res[, g]
  mean(sapply(seq_along(x), function(i) sign(cor(x[-i], y[-i]))) == sign(cor(x, y)))
}, shown$Metabolite, shown$Gene)
shown$adj_plot_support <- shown$adj_support & shown$loo_sign >= loo_sign_cut
cat("Displayed:", nrow(met_sel), "metabolites,", nrow(gene_sel), "genes,", nrow(shown), "edges (",
    sum(shown$adj_plot_support), "with TRG-adjusted support )\n")
save_table(shown %>% select(-edge_id), "Fig4D_displayed_edges")
save_table(bind_rows(met_sel %>% transmute(Node = Label, Type = "Metabolite", Class, log2FC, P = P_met),
                     gene_sel %>% transmute(Node = Gene, Type = "Gene", Class = Category, log2FC = LFC, P = DESeq2_FDR)), "Fig4D_displayed_nodes")


## ---- 4. Figure 4D --------------------------------------------------------------------------
## Three rails: pCR-enriched genes (left), metabolites (centre), non-pCR-enriched
## genes (right).  Edge colour = Spearman rho, diamond = TRG-adjusted support,
## metabolite fill = log2FC (pCR / non-pCR), gene fill = functional category.

x_left <- -1.55; x_centre <- 0; x_right <- 1.55
left  <- gene_sel %>% filter(Side == "pCR") %>% arrange(anchor_order, Category, anchor_p, desc(anchor_rho), best_p, Gene)
right <- gene_sel %>% filter(Side == "non-pCR") %>% arrange(anchor_order, Category, anchor_p, desc(anchor_rho), best_p, Gene)
y_top <- 1 + (max(nrow(left), nrow(right), nrow(met_sel)) - 1) * 0.74
left$x <- x_left; right$x <- x_right
left$y <- seq(y_top, 1, length.out = nrow(left)); right$y <- seq(y_top, 1, length.out = nrow(right))
gene_pos <- bind_rows(left, right)
met_sel$x <- x_centre; met_sel$y <- seq(y_top, 1, length.out = nrow(met_sel))

## cubic Bezier curves from metabolite to gene
curves <- shown %>%
  left_join(met_sel %>% select(Metabolite, mx = x, my = y), by = "Metabolite") %>%
  left_join(gene_pos %>% select(Gene, gx = x, gy = y), by = "Gene") %>%
  mutate(edge_order = row_number()) %>%
  uncount(51, .id = "k") %>%
  mutate(t = (k - 1) / 50, c1x = mx + 0.35 * (gx - mx), c2x = mx + 0.75 * (gx - mx),
         x = (1 - t)^3 * mx + 3 * (1 - t)^2 * t * c1x + 3 * (1 - t) * t^2 * c2x + t^3 * gx,
         y = (1 - t)^3 * my + 3 * (1 - t)^2 * t * my + 3 * (1 - t) * t^2 * gy + t^3 * gy)
fc_lim <- max(1, quantile(abs(met_sel$log2FC), 0.9))
met_sel$log2FC_plot <- pmax(pmin(met_sel$log2FC, fc_lim), -fc_lim)

p_4d <- ggplot() +
  annotate("segment", x = c(x_left, x_centre, x_right), xend = c(x_left, x_centre, x_right), y = 0.65, yend = y_top + 0.35,
           colour = "grey80", linewidth = 0.25) +
  geom_path(data = curves, aes(x, y, group = edge_id, colour = rho), linewidth = 0.35, alpha = 0.9, lineend = "round") +
  geom_point(data = filter(curves, k == 26, adj_plot_support), aes(x, y, shape = "TRG-adjusted residual-rank support"),
             size = 1.7, fill = "white", colour = "grey15", stroke = 0.4) +
  geom_point(data = met_sel, aes(x, y, fill = log2FC_plot), shape = 21, size = 4.6, colour = "grey30", stroke = 0.55) +
  geom_text(data = met_sel, aes(x, y, label = Label), nudge_y = 0.38, size = 2.9, colour = "grey15") +
  annotate("text", x = c(x_left, x_right), y = y_top + 0.85, fontface = "bold", size = 3.3,
           label = c("pCR-enriched", "non-pCR-enriched"), colour = c(response_colors[["pCR"]], response_colors[["non-pCR"]])) +
  scale_colour_gradient2(low = "#5A83AD", mid = "#EDEAE6", high = "#B86B68", limits = c(-1, 1), name = "Cor.",
                         guide = guide_colourbar(order = 3, barwidth = unit(3.5, "mm"), barheight = unit(25, "mm"))) +
  scale_fill_gradient2(low = response_colors[["non-pCR"]], mid = "white", high = response_colors[["pCR"]], limits = c(-fc_lim, fc_lim),
                       name = expression(Log[2] * " FC"), guide = guide_colourbar(order = 4, barwidth = unit(3.5, "mm"), barheight = unit(25, "mm"))) +
  scale_shape_manual(values = c("TRG-adjusted residual-rank support" = 23), name = NULL,
                     guide = guide_legend(order = 1, override.aes = list(fill = "white", colour = "grey15", size = 2.2))) +
  new_scale_fill() +
  geom_point(data = gene_pos, aes(x, y, fill = Category), shape = 21, size = 2.5, colour = "grey35", stroke = 0.35) +
  geom_text(data = filter(gene_pos, x < 0), aes(x - 0.08, y, label = Gene), hjust = 1, size = 2.9, colour = "grey20") +
  geom_text(data = filter(gene_pos, x > 0), aes(x + 0.08, y, label = Gene), hjust = 0, size = 2.9, colour = "grey20") +
  scale_fill_manual(values = category_colors, breaks = names(category_colors), drop = FALSE, name = "Category",
                    guide = guide_legend(order = 2, override.aes = list(shape = 21, size = 3.5))) +
  coord_cartesian(xlim = c(x_left - 1.1, x_right + 1.2), ylim = c(0.45, y_top + 1.2), clip = "off") +
  theme_void(base_size = 9) +
  theme(legend.position = "right", legend.title = element_text(size = 8.5), legend.text = element_text(size = 8),
        legend.spacing.y = unit(1.5, "mm"), plot.margin = margin(5, 5, 5, 5))
save_panel(p_4d, "Fig4D_metabolite_host_network", width = 7, height = max(6, 0.22 * y_top + 2.4))

message("Done: 15_Fig4D_Metabolite_host_network.R")
