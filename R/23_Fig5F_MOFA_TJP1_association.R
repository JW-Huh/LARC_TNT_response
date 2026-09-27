## =============================================================================
##  23_Fig5F_MOFA_TJP1_association.R
##  ---------------------------------------------------------------------------
##  PANEL       Figure 5F   MOFA Factor 7 score versus baseline tumour TJP1
##                          expression (20 patients with a baseline biopsy)
##
##  PURPOSE     Link the microbial Factor 7 (Figure 5B-E) to the epithelial
##              barrier gene TJP1 (ZO-1) in the matched tumour transcriptome,
##              using the full QC-passing VST matrix rather than the 2,000-gene
##              MOFA view.
##
##  INPUT       output/cache/mofa_model.RData          (script 20: mofa_input, mofa_model)
##              output/cache/mofa_associations.RData   (script 21: a = factor associations)
##              output/cache/four_omics_inputs.RData   (script 18: coherence$modalities$host$vst_qc)
##              input/metadata/patients_covariates.csv   Sex, Age, BMI (public; the
##                  restricted patients_clinical.csv is used instead when present)
##                  - without either file the partial Spearman line is omitted.
##
##  OUTPUT      output/figures/Fig5F_Factor7_TJP1  (.svg/.png)
##              output/tables/Fig5F_Factor7_TJP1_statistics.csv
##
##  METHODS     Factor 7 scores are the standardised scores of the frozen MOFA
##              model, oriented so that pCR is positive (factor_sign).  TJP1 is
##              the blind DESeq2 VST value from script 18 (all count-QC genes,
##              pooled transformation; do NOT substitute the baseline-only
##              DESeq2/VST matrix of script 10, which is fitted separately).
##              Baseline biopsies of patients in the MOFA model: n = 20.
##              Three statistics on the same patients (host_rank_screen,
##              R/_helpers_mofa.R): Pearson r, Spearman rho, and a partial
##              Spearman rho after rank regression on Age, BMI and Sex; each
##              with a permutation P value (999 permutations of the score
##              residuals, seed 20260908 + 622).  Nominal P values; a single
##              a-priori gene, so no multiplicity adjustment is applied.
##              Expected (frozen model): Pearson r = 0.48 (P = 0.029),
##              Spearman rho = 0.60 (P = 0.004), partial rho = 0.62 (P = 0.006).
##
##  WHAT CHANGES WHAT
##              gene <- "TJP1"           any other gene of the VST matrix can be
##                                       plotted with the same three tests
##              mofa_cfg$permutations    resolution of the permutation P values
##              factor                   a$axes[1]; another factor can be inspected
##
##  R PACKAGES  ggplot2
## =============================================================================

source("R/_common.R")
source("R/_helpers_mofa.R")
report_versions(c("ggplot2"))

load(cache_file("mofa_model.RData"))                                  # mofa_input, mofa_model
load(cache_file("mofa_associations.RData"))                           # a
load(cache_file("four_omics_inputs.RData"))                           # coherence

gene        <- "TJP1"                                                 # <-- host gene of interest
lead_factor <- a$axes[1]                                              # "Factor7"


## ---- 1. Baseline biopsies of the MOFA patients --------------------------------------------

host_block <- coherence$modalities$host
vst <- host_block$vst_qc                                              # biopsies x all count-QC genes (blind VST)
sm  <- host_block$sample_meta
stopifnot(identical(rownames(vst), sm$SampleID), gene %in% colnames(vst))

keep <- sm$Timepoint == "Before" & sm$SampleID %in% rownames(a$Z_all)
sm  <- sm[keep, ]; vst <- vst[keep, , drop = FALSE]
stopifnot(!anyDuplicated(sm$SubjectID))

## ORDER MATTERS for the permutation P values: the original analysis kept the
## biopsies in the column order of the count file, i.e. the alphabetical order
## of the RNA sample IDs (TNT_RNA_1, TNT_RNA_10, ...).  The correlation
## coefficients do not depend on it; the permutation P values do (slightly).
ord <- order(sm$OmicSampleID)
sm  <- sm[ord, ]; vst <- vst[ord, , drop = FALSE]
cat(sprintf("%d baseline biopsies with a %s score\n", nrow(sm), sub("Factor", "Factor ", lead_factor)))

x <- a$Z_all[match(sm$SampleID, rownames(a$Z_all)), lead_factor] * factor_sign(a, lead_factor)   # pCR direction positive
y <- vst[, gene, drop = FALSE]
patients <- read_patients()
response <- patients$Response[match(sm$SubjectID, patients$SubjectID)]


## ---- 2. Three correlation statistics on the same patients -----------------------------------

have_clinical <- has_covariates(c("Age", "BMI", "Sex"))      # patients_covariates.csv (public) or patients_clinical.csv
cv <- if (have_clinical) read_patients(clinical = TRUE)[, c("SubjectID", "Age", "BMI", "Sex")] else data.frame(SubjectID = sm$SubjectID)
cv <- cv[match(sm$SubjectID, cv$SubjectID), , drop = FALSE]

stats <- list(
  Pearson  = host_rank_screen(x, y, cv, adjust = character(), method = "pearson"),
  Spearman = host_rank_screen(x, y, cv, adjust = character(), method = "spearman")
)
if (have_clinical) {
  stats$`Partial Spearman` <- host_rank_screen(x, y, cv, adjust = c("Age", "BMI", "Sex"), method = "spearman")
} else {
  message("Age / BMI / Sex not available (input/metadata/patients_covariates.csv): the partial Spearman correlation is omitted.")
}

stats <- do.call(rbind, lapply(names(stats), function(nm) data.frame(statistic = nm, stats[[nm]][, c("n", "adjusted", "rho", "P", "t_P")])))
print(stats, row.names = FALSE, digits = 3)
save_table(stats, "Fig5F_Factor7_TJP1_statistics")


## ---- 3. Figure 5F ------------------------------------------------------------------------------

## annotation lines (plotmath so that P is italic), stacked above the legend box
symbol <- c(Pearson = "r", Spearman = "rho", `Partial Spearman` = "rho")
lines <- c(sprintf("'Baseline, n = %d'", nrow(sm)),
           sprintf("'%s %s = %.2f, '*italic(P)*' = %.3f'", stats$statistic, symbol[stats$statistic], stats$rho, stats$P))
n_lines <- length(lines)
## a text grob with explicit line positions (plotmath lines have unequal heights,
## so vjust-stacking would give uneven spacing); the block sits above the legend
annotation <- grid::textGrob(parse(text = lines), x = grid::unit(0.98, "npc"),
                             y = grid::unit(0.34, "npc") - grid::unit(seq_len(n_lines) - 1, "lines") * 1.15,
                             hjust = 1, gp = grid::gpar(fontsize = 8.5))

d <- data.frame(Score = x, Expression = y[, 1], Response = response)
p_5f <- ggplot(d, aes(Score, Expression)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey50", linewidth = 0.7) +
  geom_point(aes(colour = Response), size = 2.4) +
  annotation_custom(annotation) +
  scale_colour_manual(values = response_colors, name = NULL) +
  labs(x = paste(sub("Factor", "Factor ", lead_factor), "score"), y = sprintf("%s expression (VST)", gene)) +
  theme_tnt(base_size = 10, legend = "inside") +
  theme(aspect.ratio = 1, legend.position.inside = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = 0.3),
        legend.margin = margin(1, 5, 1, 3), legend.key.size = unit(0.4, "cm"))
save_panel(p_5f, "Fig5F_Factor7_TJP1", width = 3.8, height = 3.8)

message("Done: 23_Fig5F_MOFA_TJP1_association.R")
