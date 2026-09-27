## =============================================================================
##  16_SFig5B_Immune_module_correlations.R
##  ---------------------------------------------------------------------------
##  PANEL       Supplementary Figure 5B  fecal metabolites versus 11 tumour immune
##                                       module scores: Pearson, Spearman and
##                                       covariate-adjusted Spearman correlations
##                                       drawn as three-sector circles
##
##  PURPOSE     Check whether metabolite-immune associations are robust to the
##              choice of correlation method and to adjustment for age, sex and
##              BMI, using fixed, non-overlapping immune gene modules.
##
##  INPUT       output/cache/host_metabolite_inputs.RData   (script 14)
##              input/metadata/patients_covariates.csv (Sex, Age, BMI; public) or the
##                  restricted patients_clinical.csv - either provides the covariates
##
##  OUTPUT      output/figures/SFig5B_immune_module_correlations.{svg,png}
##              output/tables/SFig5B_immune_module_correlations.csv
##
##  METHODS     * Module score = mean of gene-wise z-scored VST expression over
##                the module's genes that pass QC (>= 2 genes required).
##              * For each metabolite (log10(conc + 1e-6)) x module: Pearson r,
##                Spearman rho, and partial Spearman = correlation of the rank
##                residuals after regressing on Age + Sex + BMI (df = n - 5).
##              * Circle radius encodes -log10 of the geometric mean of the three
##                P values - a descriptive size scale, NOT a combined test.  A
##                black outline marks a "coordinated trend" (all three
##                coefficients share the same sign).  Only SCFA, bile-acid and
##                tryptophan/indole metabolites are shown (21 in the submitted
##                panel; see the note on lithocholic acid in section 2).
##
##  R PACKAGES  dplyr, tidyr, stringr, ggplot2, ggforce, patchwork
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggforce)
  library(patchwork)
})
report_versions(c("ggplot2", "ggforce", "patchwork"))

load(cache_file("host_metabolite_inputs.RData"))                 # meta, met_log, host_expr, host_qc, immune_genes
if (!has_covariates(c("Age", "Sex", "BMI"))) {
  stop("Age, Sex and BMI (input/metadata/patients_covariates.csv or patients_clinical.csv) are required for the adjusted correlations.", call. = FALSE)
}
clinical <- read_patients(clinical = TRUE)                      # Sex / Age / BMI from whichever covariate file is present
covariates <- clinical[match(meta$SubjectID, clinical$SubjectID), c("Age", "Sex", "BMI")]
stopifnot(!anyNA(covariates))


## ---- 1. Immune modules (fixed, non-overlapping gene lists) --------------------------------

immune_modules <- list(
  NK_core        = c("KLRD1", "NCR1", "KLRK1", "KLRC1", "KLRC2", "KLRF1", "GNLY", "FGFBP2", "XCL1", "XCL2", "FCGR3A", "FCER1G", "TYROBP"),
  CD8_effector   = c("CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "RUNX3", "PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "CTSW", "CCL5", "FASLG", "TNFSF10", "IL2", "TNF", "LTA"),
  CD8_exhaustion = c("PDCD1", "TOX", "TOX2", "LAG3", "HAVCR2", "TIGIT", "ENTPD1", "LAYN", "CD244", "CD160", "NR4A1", "NR4A2", "NR4A3", "EOMES", "PRDM1"),
  Th1_core       = c("TBX21", "IFNG", "CXCR3", "IL12RB2", "IL18R1", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11"),
  Th2_core       = c("GATA3", "IL4", "IL5", "IL13", "PTGDR2", "CCR4", "IL4R", "STAT6"),
  Th17_core      = c("RORC", "IL17A", "IL17F", "CCR6", "IL23R", "KLRB1", "CCL20", "RORA", "IL21", "IL22", "IL26", "CSF2", "BHLHE40", "AHR"),
  Treg_core      = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "CCR8", "TNFRSF18", "ICOS", "LRRC32", "FGL2", "IL10", "EBI3", "NT5E"),
  Neutrophil_core = c("FCGR3B", "CEACAM8", "CSF3R", "S100A12", "FPR1", "FPR2", "CXCR1", "MMP8", "MMP9", "LTF", "CAMP", "LCN2", "OLFM4", "MNDA", "SELL"),
  MDSC_core      = c("OLR1", "ARG1", "S100A8", "S100A9", "CEBPB", "STAT3", "CXCR2", "PTGS2", "CYBB", "NCF1", "NCF2", "NCF4", "IDO1", "SLC7A2", "IL1B", "VCAN", "FCN1", "LILRB1", "CD14"),
  cDC1           = c("CLEC9A", "XCR1", "BATF3", "WDFY4", "IRF8", "CADM1", "THBD", "DNASE1L3", "SNX22", "CPVL", "CLNK"),
  TLS_Bcell      = c("CXCL13", "CCL19", "CCL21", "MS4A1", "CD79A", "CD79B", "CD74", "IGKC", "BCL6", "IL21R", "CD37", "BANK1", "MZB1", "JCHAIN")
)
stopifnot(!anyDuplicated(unlist(immune_modules)))

## module scores: mean z-score of the QC-passing member genes
scores <- sapply(immune_modules, function(g) {
  g <- intersect(g, colnames(host_expr))
  qc <- host_qc[match(g, host_qc$Gene), ]
  g <- g[ifelse(g %in% immune_genes, qc$Curated, qc$General)]
  if (length(g) < 2) return(rep(NA_real_, nrow(host_expr)))
  rowMeans(scale(host_expr[, g, drop = FALSE]))
})
module_size <- sapply(immune_modules, function(g) {
  g <- intersect(g, colnames(host_expr)); qc <- host_qc[match(g, host_qc$Gene), ]
  sum(ifelse(g %in% immune_genes, qc$Curated, qc$General))
})
scores <- scores[, colSums(is.na(scores)) == 0, drop = FALSE]
print(module_size)


## ---- 2. Three correlation coefficients per metabolite x module ------------------------------

metabolite_class <- function(x) {
  k <- tolower(gsub("_", " ", x))
  case_when(k %in% c("acetate", "propionate", "butyrate", "isobutyrate", "isovalerate", "valerate", "gamma aminobutyric acid") ~ "SCFA / related",
            str_detect(k, "cholic|deoxycholic|ursodeoxycholic|lithocholic|tauro|glyco") ~ "Bile acid",
            k %in% c("indolepropionic acid", "indole acetic acid", "indole lactic acid", "indole", "tryptamine", "tryptophan",
                     "kynurenic acid", "xanthurenic acid", "nicotinic acid") ~ "Tryptophan / indole", TRUE ~ "Other")
}

cov_df <- data.frame(Age = covariates$Age, Sex = factor(covariates$Sex), BMI = covariates$BMI)
results <- expand_grid(Module = colnames(scores), Metabolite = colnames(met_log)) %>%
  rowwise() %>%
  mutate(stats = list({
    x <- met_log[, Metabolite]; y <- scores[, Module]
    pe <- cor.test(x, y, method = "pearson"); sp <- cor.test(x, y, method = "spearman", exact = FALSE)
    rx <- residuals(lm(rank(x) ~ Age + Sex + BMI, cov_df)); ry <- residuals(lm(rank(y) ~ Age + Sex + BMI, cov_df))
    pr <- cor(rx, ry); df <- length(x) - 3 - 2
    data.frame(Pearson = unname(pe$estimate), Spearman = unname(sp$estimate), Partial = pr,
               Pearson_P = pe$p.value, Spearman_P = sp$p.value, Partial_P = 2 * pt(-abs(pr * sqrt(df / (1 - pr^2))), df), N = length(x))
  })) %>% ungroup() %>% unnest(stats) %>%
  mutate(Class = metabolite_class(Metabolite),
         Geometric_P = exp(rowMeans(log(pmax(cbind(Pearson_P, Spearman_P, Partial_P), .Machine$double.xmin)))),
         Coordinated = sign(Pearson) == sign(Spearman) & sign(Spearman) == sign(Partial) & Pearson != 0)
save_table(results, "SFig5B_immune_module_correlations")

classes <- c("SCFA / related", "Bile acid", "Tryptophan / indole")               # <-- metabolite classes shown

## NOTE  In the laboratory table used for the submission the lithocholic acid
## column was spelled "Lithocholi_acid", so the class regex above did not
## recognise it as a bile acid and the submitted panel shows 21 metabolites.
## The curated input corrects the name; to reproduce the submitted panel the
## metabolite is excluded here explicitly.  Set `excluded_metabolites <-
## character()` to draw all 22 class members (a 9th bile-acid row appears).
excluded_metabolites <- c("Lithocholic_acid")
results <- results %>% filter(Class %in% classes, !Metabolite %in% excluded_metabolites)

## metabolite order: Ward clustering of the coefficient profiles within each class
order_metabolites <- unlist(lapply(classes, function(cl) {
  wide <- results %>% filter(Class == cl) %>% select(Metabolite, Module, Pearson, Spearman, Partial) %>%
    pivot_wider(names_from = Module, values_from = c(Pearson, Spearman, Partial))
  x <- as.matrix(wide[, -1]); x[!is.finite(x)] <- 0
  wide$Metabolite[if (nrow(x) > 1) hclust(dist(x), method = "ward.D2")$order else 1]
}))


## ---- 3. Supplementary Figure 5B ----------------------------------------------------------
## Each circle has three sectors (clockwise from 12 o'clock): partial Spearman,
## Pearson, Spearman.  Radius = -log10(geometric-mean P) capped at 3.

results <- results %>%
  mutate(x = match(Module, colnames(scores)), y = match(Metabolite, rev(order_metabolites)),
         radius = 0.10 + 0.30 * sqrt(pmin(-log10(Geometric_P), 3) / 3),
         Outline = ifelse(Coordinated, "Coordinated trend", "Single/mixed trend"))
sectors <- results %>%
  select(Module, Metabolite, Class, x, y, radius, Outline, Partial, Pearson, Spearman, Partial_P, Pearson_P, Spearman_P) %>%
  pivot_longer(c(Partial, Pearson, Spearman), names_to = "Method", values_to = "Correlation") %>%
  mutate(P = case_when(Method == "Partial" ~ Partial_P, Method == "Pearson" ~ Pearson_P, TRUE ~ Spearman_P),
         start = c(Partial = 0, Pearson = 2 * pi / 3, Spearman = 4 * pi / 3)[Method], end = start + 2 * pi / 3,
         Symbol = case_when(P < 0.001 ~ "$", P < 0.01 ~ "#", P < 0.05 ~ "*", TRUE ~ ""),
         sx = x + 0.58 * radius * sin((start + end) / 2), sy = y + 0.58 * radius * cos((start + end) / 2))
module_labels <- paste0(gsub("_", " ", colnames(scores)), " (", module_size[colnames(scores)], ")")
class_breaks <- results %>% distinct(Metabolite, Class) %>% mutate(y = match(Metabolite, rev(order_metabolites))) %>%
  group_by(Class) %>% summarise(ymin = min(y) - 0.5, ymax = max(y) + 0.5, .groups = "drop")

p_circles <- ggplot(sectors) +
  geom_rect(data = class_breaks, aes(xmin = 0.4, xmax = ncol(scores) + 0.6, ymin = ymin, ymax = ymax), fill = NA, colour = "grey40", linewidth = 0.4) +
  geom_text(data = class_breaks, aes(x = ncol(scores) + 0.9, y = (ymin + ymax) / 2, label = Class), angle = -90, size = 3, fontface = "bold") +
  geom_arc_bar(aes(x0 = x, y0 = y, r0 = 0, r = radius, start = start, end = end, fill = Correlation, colour = Outline), linewidth = 0.4) +
  geom_text(aes(sx, sy, label = Symbol), size = 2) +
  scale_x_continuous(breaks = seq_len(ncol(scores)), labels = module_labels) +
  scale_y_continuous(breaks = seq_along(order_metabolites), labels = gsub("_", " ", rev(order_metabolites))) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "Correlation") +
  scale_colour_manual(values = c("Coordinated trend" = "black", "Single/mixed trend" = "grey65"), name = NULL) +
  coord_fixed(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_tnt(base_size = 9, legend = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank(),
        plot.margin = margin(5, 30, 5, 5))

## key: sector layout, size scale and symbols
key <- data.frame(start = c(0, 2 * pi / 3, 4 * pi / 3), end = c(2 * pi / 3, 4 * pi / 3, 2 * pi))
size_key <- data.frame(y = c(-2, -3.2, -4.4), value = 1:3, radius = 0.10 + 0.30 * sqrt((1:3) / 3))
p_key <- ggplot() +
  geom_arc_bar(data = key, aes(x0 = 0, y0 = 0, r0 = 0, r = 0.5, start = start, end = end), fill = "white", colour = "grey40") +
  annotate("text", x = c(0.6, 0, -0.6), y = c(0.35, -0.8, 0.35), hjust = c(0, 0.5, 1),
           label = c("Partial Spearman", "Pearson", "Spearman"), size = 2.8) +
  geom_circle(data = size_key, aes(x0 = 0, y0 = y, r = radius), fill = "white", colour = "grey30") +
  geom_text(data = size_key, aes(x = 0.8, y = y, label = value), size = 2.8) +
  annotate("text", x = 0, y = -1.3, label = '-Log[10]~"(geometric mean"~italic(p)*")"', parse = TRUE, size = 2.8) +
  annotate("text", x = -1.2, y = c(-5.2, -5.6, -6.0), hjust = 0, size = 2.8, parse = TRUE,
           label = c('"*"~italic(p) < 0.05', '"#"~italic(p) < 0.01', '"$"~italic(p) < 0.001')) +
  coord_fixed(xlim = c(-2.7, 2.7), ylim = c(-6.2, 1.2), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(5, 25, 5, 5))

save_panel(p_circles + p_key + plot_layout(widths = c(1, 0.5)), "SFig5B_immune_module_correlations", width = 8.8, height = 7)
cat("Displayed:", length(order_metabolites), "metabolites x", ncol(scores), "modules\n")

message("Done: 16_SFig5B_Immune_module_correlations.R")
