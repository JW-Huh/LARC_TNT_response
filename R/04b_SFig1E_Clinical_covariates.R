## =============================================================================
##  04b_SFig1E_Clinical_covariates.R
##  ---------------------------------------------------------------------------
##  PANEL       Supplementary Figure 1E  variance in baseline microbiome
##              composition explained by each clinical covariate (PERMANOVA)
##
##  PURPOSE     Put the response association (S1D, R^2 = 4.4 %) into context:
##              how much of the between-patient Bray-Curtis variation is
##              explained by each of 26 clinical variables tested one at a time?
##
##  INPUT       input/microbiome/metaphlan_sgb.csv
##              input/metadata/stool_samples.csv, patients.csv
##              input/metadata/patients_clinical.csv   *** RESTRICTED ***
##                  patient-level clinical covariates (26 variables); available
##                  from the authors on request, not part of the public deposit.
##              input/metadata/sfig1E_permanova_frozen.csv
##                  the per-covariate PERMANOVA results (R2, F, P, N) of the
##                  submitted run - no patient-level information.
##
##  TWO MODES   with the clinical table   every PERMANOVA is recomputed and
##                  compared with the frozen results.
##              without it (public copy)  the panel is drawn from the frozen
##                  results, so the figure is still reproduced.
##
##  OUTPUT      output/figures/SFig1E_clinical_PERMANOVA.{svg,png}
##              output/tables/SFig1E_clinical_PERMANOVA.csv
##
##  METHODS     For every covariate: keep the patients with a recorded value,
##              compute Bray-Curtis dissimilarity on their SGB profiles and fit
##              a one-variable PERMANOVA (vegan 2.6-6.1, adonis2(), 999
##              permutations, seed 123).  Character covariates (X/O,
##              Past/Present/X, M/F) are categorical; numeric covariates enter
##              as a single linear term.  R^2 (%) is displayed.  These are 26
##              separate marginal tests, not a jointly adjusted model, and the
##              denominators (N) differ because of missing values (20-26;
##              empty cells of the clinical table are read as NA).
##
##  R PACKAGES  vegan, ggplot2
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages(library(vegan))
report_versions(c("vegan", "ggplot2"))


## ---- 1. Clinical covariates ---------------------------------------------------

frozen <- read_csv_na(input_file("metadata", "sfig1E_permanova_frozen.csv"))   # results of the submitted run
have_clinical <- file.exists(input_file("metadata", "patients_clinical.csv"))
if (!have_clinical) message("patients_clinical.csv (restricted) not found: Supplementary Figure 1E is drawn from the frozen PERMANOVA results.")

## Variables tested, grouped as in the panel.  <-- add / remove covariates here
covariate_groups <- list(
  "Patient-related"  = c("Sex", "Age", "BMI", "ASA", "Smoking", "Drinking", "Cancer_family_history"),
  "Tumor-related"    = c("TRG_score", "cT_stage", "cN_stage", "CEA"),
  "Comorbidities"    = c("Diabetes", "Hypertension", "Heart_disease", "Pulmonary_disease",
                         "Liver_disease", "Cerebrovascular_accident"),
  "Blood biomarkers" = c("Glucose", "Albumin", "Hemoglobin", "CRP", "Segmented_neutrophil",
                         "Lymphocyte", "NLR", "WBC", "Creatinine")
)
covariates <- unlist(covariate_groups, use.names = FALSE)

## Axis labels for the panel (the submitted panel spells "Creatinine" as "Creatine")
covariate_labels <- setNames(gsub("_", " ", covariates), covariates)
covariate_labels[c("cT_stage", "cN_stage")] <- c("cT stage", "cN stage")


## ---- 2. One-variable PERMANOVA per covariate -----------------------------------

permanova_one <- function(field) {
  keep <- !is.na(baseline[[field]])
  value <- baseline[[field]][keep]
  stopifnot(sum(keep) >= 4, length(unique(value)) > 1)
  d <- vegdist(sgb[keep, , drop = FALSE], method = "bray")     # complete cases only
  set.seed(123)                                                # <-- permutation seed
  fit <- adonis2(d ~ value, permutations = 999)
  data.frame(Covariate = field, N = sum(keep),
             Type = ifelse(is.numeric(value), "numeric", "categorical"),
             R2_percent = 100 * fit$R2[1], F = fit$F[1], P = fit$`Pr(>F)`[1])
}

if (have_clinical) {
  baseline <- read_stool_samples(timepoint = "Before", clinical = TRUE)
  sgb <- read_metaphlan("sgb", samples = baseline$SampleID, drop_empty = FALSE)   # 26 x SGBs (%)
  stopifnot(all(covariates %in% names(baseline)))
  results <- do.call(rbind, lapply(covariates, permanova_one))
  results$Group <- rep(names(covariate_groups), lengths(covariate_groups))
  ## check against the frozen results of the submitted run
  chk <- merge(results, frozen, by = "Covariate", suffixes = c("", "_frozen"))
  cat(sprintf("Recomputed vs frozen: max |R2 diff| = %.2e %%, max |P diff| = %.3f\n",
              max(abs(chk$R2_percent - chk$R2_percent_frozen)), max(abs(chk$P - chk$P_frozen))))
} else {
  results <- frozen[frozen$Covariate %in% covariates, ]
}
results <- results[order(-results$R2_percent), ]
print(results, row.names = FALSE, digits = 3)
save_table(results, "SFig1E_clinical_PERMANOVA")


## ---- 3. Supplementary Figure 1E -------------------------------------------------
## Design (as submitted): vertical bars in decreasing R^2, "R^2 = x.xx" written
## vertically inside each bar in white, covariate names rotated 90 degrees,
## group legend on the right.

group_colors <- c("Blood biomarkers" = "#C98845", "Comorbidities" = "#A1669E",
                  "Patient-related" = "#5FAD93", "Tumor-related" = "#6985BE")

results$Covariate <- factor(results$Covariate, levels = results$Covariate)

p_s1e <- ggplot(results, aes(Covariate, R2_percent, fill = Group)) +
  geom_col(width = 0.8) +
  geom_text(aes(y = R2_percent / 2, label = sprintf("R² = %.2f", R2_percent)),
            angle = 90, colour = "white", size = 2.9) +
  scale_fill_manual(values = group_colors) +
  scale_x_discrete(labels = covariate_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = seq(0, 9, 3)) +
  labs(x = NULL, y = "Variance explained (%)", fill = "Group") +
  theme_tnt(base_size = 10, legend = "right") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
save_panel(p_s1e, "SFig1E_clinical_PERMANOVA", width = 9, height = 3.4)

message("Done: 04b_SFig1E_Clinical_covariates.R")
