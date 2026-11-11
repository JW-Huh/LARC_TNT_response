# LARC_TNT_response — reproducible analysis code

Gut microbiome, fecal metabolome and tumour transcriptome of locally advanced
rectal cancer (LARC) patients treated with total neoadjuvant therapy (TNT):
pathological complete response (pCR, TRG 0, n = 11) versus non-pCR (TRG 1–3,
n = 15).

This folder is a self-contained R project that regenerates **every panel of
Figures 1–5 and Supplementary Figures 1–6** of the manuscript, the numbers
quoted in the text, and the tables behind each panel, from curated input
files in `input/`.

```
LARC_TNT_response/
├── LARC_TNT_response.Rproj   open this in RStudio (sets the working directory)
├── README.md                 this file
├── run_all.R                 runs every script in order
├── R/                        28 analysis scripts + 3 shared helper files
├── input/                    curated, ready-to-use input files (see input/README.md)
│   ├── metadata/  microbiome/  metabolome/  host_rnaseq/  mofa/  random_forest/
│   └── INPUT_MANIFEST.csv    file list with row / column counts and MD5 checksums
├── output/                   created by the scripts (not distributed)
│   ├── figures/              one .svg (editable) + one .png (quick look) per panel
│   ├── tables/               the numbers behind each panel (.csv)
│   └── cache/                intermediate R objects passed between scripts
└── dev/                      prepare_inputs_from_source.R: how input/ was exported
                              from the raw laboratory tables (provenance only)
```

---

## 1. Quick start

```r
# 1. open LARC_TNT_response.Rproj in RStudio  (or setwd() to this folder)
# 2. install the packages listed in section 3
# 3. run everything (about 15 min without refits):
source("run_all.R")
# ... or run a single script:
source("R/01_Fig1BC_SFig1ABD_Diversity.R")
```

Every script starts with `source("R/_common.R")`, which checks that R is
running in the package root, creates `output/`, and defines the colour
palettes, the ggplot2 theme and a few helper functions used everywhere.
Scripts that depend on the output of an earlier script say so in their
header (`INPUT  output/cache/...`); `run_all.R` runs them in the right order.

Figures are written as **SVG** (vector, for the final figure assembly in
Illustrator/Inkscape) and **PNG** (300 dpi, for quick inspection). The
submitted figures were assembled from the SVG panels; panel letters, the
final font and some manual spacing were added during assembly and are not
part of the R output.

### Restricted inputs

Two files contain patient-level clinical or identifying data and are **not in
the public repository** (`restricted = TRUE` in `input/INPUT_MANIFEST.csv`):

| file | content | used by | without it |
|---|---|---|---|
| `input/metadata/patients_clinical.csv` | 26 clinical variables per patient (stage, comorbidities, laboratory values, …) | 04b | the public file `input/metadata/patients_covariates.csv` (Sex, Age, BMI) provides the only covariates used by the adjusted correlations (16, 22, 23), and script 04b draws Supplementary Figure 1E from the frozen per-covariate PERMANOVA results (`sfig1E_permanova_frozen.csv`). The full table is available from the corresponding author on reasonable request. |
| `input/metadata/id_crosswalk_restricted.csv` | study ID ↔ laboratory / hospital identifiers | `dev/prepare_inputs_from_source.R` only | not needed for any analysis |

Everything else — including the de-identified gene-level count matrix of the
27 tumour biopsies (`input/host_rnaseq/tumor_rna_counts.csv`, the same matrix
deposited in GEO GSE339362) — is in the repository, so `run_all.R` reproduces
every figure. If the count matrix is removed, `run_all.R` skips scripts 10 and
13–23 (Figures 4C–E, 5A–F, S5, S6) and everything else still runs.

---

## 2. Script index

Scripts are numbered in execution order. "Frozen" means the script replays a
stored result of the submitted analysis instead of recomputing a stochastic or
very slow step; each of those scripts has a `refit_* <- TRUE` switch at the top
that turns the full computation back on (see section 5).

| script | panels | what it does | needs |
|---|---|---|---|
| `_common.R` | – | paths, colours, theme, readers (`read_patients`, `read_stool_samples`, `read_metaphlan`, `read_metabolites`, …), taxon-name formatting, Wilcoxon / Hedges' g helpers, figure export | – |
| `_helpers_enrichment_map.R` | – | keyword rules that name the enrichment-map clusters (Figure 4A refit only) | – |
| `_helpers_mofa.R` | – | MOFA input construction, mixed-model / permutation statistics, feature selection, partial correlations (scripts 20–23) | – |
| `00_Fig1A_Fig3E_SFig4A_Study_schematics.R` | 1A, 3E, S4A | data-availability UpSet plot (data-driven) and the three schematic diagrams | metadata |
| `01_Fig1BC_SFig1ABD_Diversity.R` | 1B, 1C, S1A, S1B, S1D | alpha diversity (evenness, inverse Simpson, Shannon, richness), Bray–Curtis PCoA + PERMANOVA | `metaphlan_sgb.csv` |
| `02_Fig1D_SFig1C_Taxonomic_composition.R` | 1D, S1C | top-10 family composition, by group and by sample | `metaphlan_family.csv` |
| `03_Fig1EFG_SFig1FGH_Species_ecology.R` | 1E, 1F, 1G, S1F, S1G, S1H | microbial-group radar, response-associated species dumbbell + annotation tiles, co-abundance network, species heatmaps, CRC-SGB change after RT | `metaphlan_species.csv`, `metaphlan_sgb.csv`, `microbial_groups.csv` |
| `04_Fig1H_SGB_PLS_DA.R` | 1H | PLS-DA of SGB profiles, VIP scores | `metaphlan_sgb.csv` |
| `04b_SFig1E_Clinical_covariates.R` | S1E | one-variable PERMANOVA of 26 clinical covariates (recomputed with the restricted table, otherwise replayed from the frozen results) | `sfig1E_permanova_frozen.csv` |
| `05_Fig2ABC_KO_categories.R` | 2A, 2B, 2C | KO-level Wilcoxon tests, category perturbation score vs random KO sets, KO heatmap | `humann_ko_cpm.tsv`, `ko_functional_categories.csv` |
| `06_Fig2DEF_CAZyme_analysis.R` | 2D, 2E, 2F | CAZyme class shift, subfamily volcano, starch-capacity score | `cazy_*_tpm.tsv` |
| `06b_SFig2_MetaCyc_pathways.R` | S2A–M | MetaCyc pathway heatmap, sucrose/rhamnose balance groups vs TRG, genus / species contributions | `humann_pathway_*.tsv` |
| `07_Fig3BF_SFig3A_Metabolite_effects.R` | 3B, 3F, S3A | Hedges' g forest plots with bootstrap intervals, tryptophan-branch ternary display | `fecal_metabolites.csv` |
| `08_Fig3ACD_SFig3B_Metabolite_associations.R` | 3A, 3C, 3D, S3B | metabolome PCoA + envfit vectors, species–SCFA correlations, *E. rectale* EC 2.3.1.9 vs butyrate, bile-acid indices | metabolome + `metaphlan_species.csv`, `humann_ec_2.3.1.9_species_stratified.tsv` |
| `09a_SFig4_Random_forest_inputs.R` | – | loads / checks the eight frozen assay matrices for the random forests | `rf_inputs_frozen.rds` |
| `09b_SFig4_Random_forest_models.R` | (S4B–J) | **frozen**: replays and re-verifies the submitted cross-validated forests; `refit_rf <- TRUE` retrains (1–2 h) | `rf_results_submitted.rds` |
| `09c_SFig4_Random_forest_ROC_SHAP.R` | S4B–J | ROC curves and SHAP / importance panels | cache of 09b |
| `10_Fig4_SFig5A_Host_RNAseq_preparation.R` | S5A | DESeq2 (pCR vs non-pCR) and VST of the 20 baseline biopsies, volcano plot | `tumor_rna_counts.csv` |
| `11_Fig4B_Host_gene_set_enrichment.R` | 4B | **frozen**: 22 selected gene sets from the stored fgsea table (`refit_gsea <- TRUE` reruns fgsea) | `gsea_targeted_results_frozen.csv` |
| `12_Fig4A_Host_enrichment_map.R` | 4A | **frozen**: enrichment map re-screened from the stored global GSEA table, clusters from the frozen node table | `gsea_global_results_frozen.csv`, `enrichment_map_nodes_frozen.csv` |
| `13_Fig4C_Species_host_correlations.R` | 4C | species × 18 host genes Spearman heatmap with gene-set membership | cache of 10 + 11, `metaphlan_species.csv` |
| `14_Fig4DE_SFig5B_Cross_omics_preparation.R` | – | metabolite / host-gene QC for the 17 patients with both layers | cache of 10 |
| `15_Fig4D_Metabolite_host_network.R` | 4D | metabolite–gene rail network with TRG-adjusted support | cache of 14 |
| `16_SFig5B_Immune_module_correlations.R` | S5B | immune-module scores vs metabolites: Pearson / Spearman / partial Spearman | cache of 14, `patients_covariates.csv` |
| `17_Fig4E_Species_metabolite_host_triads.R` | 4E | six species–metabolite–gene triads | cache of 14, `legacy_paired_subjects.csv` |
| `18_Fig5_Four_omics_preparation.R` | – | matched species / KO / metabolite / host blocks, PCoA scores for all six layer pairs | metagenome + metabolome + `tumor_rna_counts.csv` |
| `19_Fig5A_Procrustes_GPA.R` | 5A | pairwise Procrustes and generalized Procrustes analysis, both tested with the same patient-block permutation scheme (unrestricted PROTEST P kept as `P_unrestricted`) | cache of 18 |
| `20_Fig5_MOFA_preparation_and_training.R` | – | rebuilds the MOFA views, verifies them against the frozen training input, loads the frozen model (`refit_mofa <- TRUE` retrains with MOFA2) | cache of 18, `mofa_*_frozen.rds` |
| `21_Fig5BCD_MOFA_response_associations.R` | 5B, 5C, 5D | factor–response associations (mixed model, Hedges' g, patient permutations), variance explained, Factor 7 × Factor 11 plane | cache of 20 |
| `22_Fig5E_SFig6_MOFA_feature_network.R` | 5E, S6 | Factor 7 feature loadings; adjusted partial-correlation network of those features | cache of 20 + 21, `patients_covariates.csv` (S6) |
| `23_Fig5F_MOFA_TJP1_association.R` | 5F | Factor 7 score vs baseline *TJP1* expression (Pearson, Spearman, partial Spearman) | cache of 18 + 20 + 21 |

Each script header repeats this information in more detail: **PANELS**,
**PURPOSE**, **INPUT**, **OUTPUT**, **METHODS** (algorithm, package version and
a short description of the principle), and, where useful, a **WHAT CHANGES
WHAT** list that tells the reader which parameter changes which part of the
figure. Tunable parameters are marked in the code with `# <--`.

---

## 3. Software

The verification runs used **R 4.3.2** (Windows, UCRT) with the package
versions below. Newer versions should work; the versions matter only when
you want to reproduce permutation P values to the last digit. The re-run
with the corrected metabolite table (November 2026, section 7) used
**R 4.6.1** with the package versions current at that time; the result
tables of all analyses that do not involve the metabolome agreed with the
R 4.3.2 results (largest relative difference 1 × 10⁻⁵, in the DESeq2
statistics; see section 7 for what differs in the drawings).

*Wilcoxon tests and the R version.* From R 4.6.0 `stats::wilcox.test`
returns an exact conditional P value for tied data by default, whereas
earlier versions use the normal approximation whenever there are ties.
`wilcox_p()` in `R/_common.R` therefore states the rule of the submitted
analysis explicitly (exact P without ties and n < 50 per group, otherwise
normal approximation with continuity correction) and returns P with 12
significant digits, so that every R version gives the same P values and the
same order of tied P values.

| package | version | used for |
|---|---|---|
| ggplot2 | 4.0.3 | all figures |
| patchwork | 1.3.2 | panel composition, insets |
| dplyr / tidyr / stringr | 1.1.4 / 1.3.x / 1.5.x | data wrangling |
| ggbeeswarm, ggrepel, ggforce, ggnewscale, svglite | 0.7.2, 0.9.6, 0.5.0, 0.5.2, ≥ 2.1 | plot details, SVG export |
| vegan | 2.6-6.1 | diversity, Bray–Curtis, PERMANOVA, envfit, Procrustes |
| ape | ≥ 5.7 | PCoA (Figure 5) |
| mixOmics | 6.26.0 | PLS-DA (Figure 1H) |
| ComplexHeatmap / circlize | 2.18.0 / 0.4.18 | heatmaps (1F tiles, 2C, S1F/G, S2A, 4C) |
| igraph / ggraph | 2.3.2 / 2.2.2 | networks (1G, 4A, S6) |
| pROC | 1.19.0.1 | ROC curves |
| randomForest, rsample, fastshap | 4.7-1.2, 1.3.0, 0.1.1 | random-forest refit only |
| DESeq2 (Bioconductor) | 1.42.1 | host RNA-seq (10, 18) |
| fgsea (Bioconductor), msigdbr | 1.28.0, MSigDB 2026.1.Hs | GSEA refit only |
| nlme | 3.1-x | mixed models (MOFA associations) |
| MOFA2 (Bioconductor) + mofapy2 via reticulate | 1.12.1 | MOFA refit only |

```r
install.packages(c("ggplot2", "patchwork", "dplyr", "tidyr", "stringr", "ggbeeswarm",
                   "ggrepel", "ggforce", "ggnewscale", "svglite", "vegan", "ape",
                   "igraph", "ggraph", "pROC", "nlme", "randomForest", "rsample", "fastshap"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("mixOmics", "ComplexHeatmap", "DESeq2", "fgsea"))   # + "MOFA2" for a MOFA refit
```

---

## 4. Inputs

`input/README.md` documents every file and column. In short:

* **metadata/** — `patients.csv` (26 patients: response, TRG), `stool_samples.csv`
  (42 metagenomes → patient, visit), `rna_samples.csv` (27 biopsies),
  `patients_covariates.csv` (Sex, Age, BMI), `sfig1E_permanova_frozen.csv`,
  `legacy_paired_subjects.csv`; the full clinical table is restricted.
* **microbiome/** — MetaPhlAn 4 profiles at five ranks, HUMAnN 3 KO / EC /
  MetaCyc tables, dbCAN CAZyme TPM, the KO → functional-category table
  (Supplementary Table 3) and the curated microbial-group memberships
  (Supplementary Table 1).
* **metabolome/** — 40 targeted fecal metabolites (µmol/g) for 32 samples
  (corrected laboratory report of November 2026, non-detected = 0; see
  section 7). The laboratory result files are in MetaboLights (MTBLS16367:
  bile acids and tryptophan-related compounds; MTBLS16378: SCFAs).
* **host_rnaseq/** — the gene-level count matrix of the tumour biopsies (as
  deposited in GEO), frozen GSEA results and the enrichment-map node table,
  and the MSigDB gene sets used.
* **random_forest/**, **mofa/** — frozen inputs and results of the two
  stochastic / slow analyses (see section 5).

Sample identifiers are anonymised (`S22207`, `Sample_7`, `TNT_RNA_44`); the
crosswalk to laboratory identifiers is restricted.

**Row order matters for permutation tests.** Permutation and bootstrap
procedures draw random numbers per row, so the same seed gives slightly
different P values if the samples are in a different order. The columns
`AnalysisOrder` (stool_samples.csv, fecal_metabolites.csv) and
`EnrolmentOrder` (patients.csv) reproduce the sample order of the original
analysis, and the readers in `_common.R` return rows in that order. If you
build your own sample tables, keep these columns.

**supplementary_tables/** holds Supplementary Table 2 of the manuscript as
submitted (`Supplementary Table 2.xlsx`, sheets S2A–S2H: response-group
comparisons of SGBs, species, genera, ECs, KOs, MetaCyc pathways,
metabolites and tumour genes at baseline and after RT). It is provided here
because the workbook is too large to be merged into the reviewer PDF; the
same statistics are recomputed by scripts 03, 05, 06b, 07 and 10 (see the
note on the Wilcoxon variant in section 6).

---

## 5. Frozen results and refit switches

Four analyses are either stochastic without a reproducible random state or
take hours; for these the submitted result is shipped in `input/` and the
script replays and *re-verifies* it by default:

| analysis | frozen file | replay checks | switch |
|---|---|---|---|
| random forests (S4) | `random_forest/rf_results_submitted.rds` | AUC of every model recomputed from the stored out-of-fold predictions | `refit_rf <- TRUE` in 09b (1–2 h; SHAP draws differ → a new fit) |
| targeted GSEA (4B) | `host_rnaseq/gsea_targeted_results_frozen.csv` | the 22 displayed NES / P / FDR are compared with the submitted values (tolerance 1e-8) | `refit_gsea <- TRUE` in 11 |
| global GSEA + enrichment map (4A) | `host_rnaseq/gsea_global_results_frozen.csv`, `enrichment_map_nodes_frozen.csv` | node selection and edges are recomputed from the table; only cluster membership / labels are taken from the frozen nodes | `refit_gsea <- TRUE` in 12 (needs msigdbr) |
| MOFA (5B–F, S6) | `mofa/mofa_input_frozen.rds`, `mofa_trained_model_frozen.rds` (retrained in November 2026 on the corrected metabolite table, section 7) | the four views are rebuilt from the raw data and must equal the frozen training input (tolerance 1e-7); variance explained is reconstructed from scores × weights and compared with the stored values | `refit_mofa <- TRUE` in 20 (needs MOFA2 + Python; `freeze_refit <- TRUE` stores the refit as the new frozen model) |

All other analyses are recomputed from the raw tables every time.

---

## 6. Verification against the submitted figures

All 28 scripts were run from a clean `output/` folder (R 4.3.2, Windows).
Every numeric result printed in a panel or quoted in the text was checked.
The list below is the record of the submitted analysis (first laboratory
report of the metabolome); the values that changed with the corrected
metabolite table are listed in section 7:

* diversity P values (evenness 0.047, inverse Simpson 0.087, Shannon 0.148,
  richness 0.775), PERMANOVA S1D (R² 0.044, P 0.291), S1E R² values;
* PLS-DA VIPs (1H), KO category Z / P (2B), KO stars (2C), CAZyme P values
  (2D–F), MetaCyc pathway P values and group tests (S2);
* Hedges' g and bootstrap intervals of all seven metabolites of 3B, ratio
  test P = 0.023 (3F), metabolome PERMANOVA (R² 0.134, P 0.063), *E. rectale*
  EC 2.3.1.9 rho 0.52 / P 0.020 (3D), bile-acid index P values (S3B);
* random-forest AUCs (S4B), the 22 gene-set statistics (4B), Fig 4C stars,
  the 4D / 4E node and edge sets;
* pairwise Procrustes r 0.715 / 0.487 / 0.440 with their intervals and
  unrestricted PROTEST P values (0.0001 / 0.0015 / 0.1403, column
  `P_unrestricted` in `Fig5A_pairwise_procrustes.csv`), GPA agreement 0.774,
  mean pairwise r 0.502, P 0.0007 (5A); MOFA
  Factor 7 P 0.011 (BH q 0.154, max-F P 0.163), Factor 11 P 0.083, plane
  R² 0.157 / P 0.020 (5B–D); joint R² 0.661; the 12 Factor 7 features and
  their loadings (5E); TJP1 Pearson r 0.48 (P 0.028), Spearman 0.60
  (P 0.006), partial 0.62 (P 0.006) (5F).

Permutation P values are reproduced to the last printed digit wherever the
sample order of the original run is known (see `AnalysisOrder` /
`EnrolmentOrder`, section 4).

Per-feature Wilcoxon P values (2A/2C, S2) are computed with
`stats::wilcox.test` (exact when there are no ties, otherwise the normal
approximation with continuity correction). Supplementary Table 2 was
generated with the exact conditional test, which handles ties through the
permutation distribution of the rank sum (`coin::wilcox_test(distribution =
"exact")`); the two agree exactly for features without ties and differ
slightly for tied features (e.g. K02406 flagellin, baseline: 0.0069 here vs
0.0054 in Table S2E). Feature selections and significance stars are the same
under both variants.

Changes made during revision (the code now differs from the submitted PDF
on purpose):

| panel | change | where |
|---|---|---|
| 2C | significance stars use the thresholds of the submitted panel and legend (`***` P < 0.005, `**` P < 0.01); the first public release used 0.001 / 0.01, which showed K13542 and K20491 as `**` | 05, `star_cuts` |
| 5A insets | the pairwise Procrustes P values are from the patient-block permutation test (patients exchanged within identical visit-availability patterns, 9,999 permutations), the same null as the GPA test, instead of unrestricted PROTEST permutations; r, intervals and the GPA statistics are unchanged, and the unrestricted P is kept in the table | header of 19 |

Known, documented differences from the submitted PDF panels:

| panel | difference | where explained |
|---|---|---|
| 1A, 3E, S4A | the timeline, the tryptophan-pathway scheme and the workflow scheme are drawn with plain boxes; the submitted versions were finished in Illustrator | header of 00 |
| S1C | the submitted legend lists three families (Enterobacteriaceae, Lactobacillaceae, Selenomonadaceae) that are not in the baseline top-10 (legend artefact of the after-RT panel); reproduced optionally | 02, `s1c_legend_extra` |
| S2A | the colour bar is labelled log10 ratio (the submitted legend said Log2FC; the values are log10) | header of 06b |
| 1G, 4A, S6 | network node positions come from a seeded force-directed layout; the submitted panels were arranged by hand (1G reproduces the Cytoscape coordinates where available) | 03, 12, 22 |
| S5B | lithocholic acid is excluded explicitly: in the laboratory table used for the submission its column was misspelled (`Lithocholi_acid`) and therefore not recognised as a bile acid; the curated input corrects the name | 16, `excluded_metabolites` |
| 4B | the first gene set ("PI3K–AKT signaling with NF-κB") is the KEGG MEDICUS set N01339 (NNK/NNN → PI3K signaling); the display name is the submitted one | 11, `catalog` |
| 5B–F, S6 | the rebuilt MOFA views are compared with the frozen training input (relative tolerance 10⁻⁴; in the submitted model the species, metabolite and host views were identical and the KO view agreed within 1.5 × 10⁻⁵ relative, the same 2,000 features); the frozen model is used, because a retrain (`refit_mofa`) can change factor numbering and sign | 20 |
| 5E | the lollipop shows the normalised loading (sign oriented towards pCR) without the Hedges' g colour of the original draft script; g is in the table | 22 |

Design: colours, shapes, axis labels and panel proportions follow the
submitted PDFs (colour codes were read from the PDFs). Fonts differ slightly
(the PDFs use Arial; R uses the system sans-serif) and a few labels in dense
panels are placed by `ggrepel` and may land a few millimetres from where they
were in the submission.

---

## 7. Corrected metabolite table (November 2026)

In November 2026 the laboratory re-issued the targeted metabolomics report.
In the first report, cells of analytes that had not been detected in a
sample had been filled in with values. 171 of the 1,280 cells (32 samples ×
40 analytes) changed: 131 values became non-detects (hyodeoxycholic acid 31,
ursodeoxycholic acid 29, tauroursodeoxycholic/taurohyodeoxycholic acid 28,
taurochenodeoxycholic acid 18, indole-3-lactic acid 8, histamine 5,
indole-3-acetic acid 4, and one or two cells of six other analytes) and 40
values were re-quantified (histamine 27, taurodeoxycholic acid 12,
tryptamine 1). `input/metabolome/fecal_metabolites.csv` holds the corrected
values (non-detected = 0); `dev/update_fecal_metabolites_v2.py` shows how
the file was rebuilt, and Supplementary Table 2 (sheet S2G) was recomputed.

All scripts were re-run with the corrected table (R 4.6.1). Results that do
not involve the metabolome did not change. Results that changed:

| panel | result | first report | corrected table |
|---|---|---|---|
| 3A | PCoA axes; PERMANOVA | 64.6 % / 15.6 %; R² 0.134, P 0.063 | 64.5 % / 15.7 %; R² 0.134, P 0.063 |
| 3B | metabolites shown | 7 | 6: indole-3-lactic acid is no longer different between the groups (P 0.038 → 0.385; g 0.40 [−0.48, 1.08]); g and intervals of the other six are unchanged |
| S3A | metabolites tested | 35 | 34 (histamine added; ursodeoxycholic and tauroursodeoxycholic/taurohyodeoxycholic acid no longer testable) |
| 3F | reductive / oxidative branch index, P | 0.023 | 0.033 |
| S3B | bile-acid indices, P | secondary/primary 0.15; CA/CDCA-derived 0.60; total 7α 0.26; CDCA-derived 7α 0.37 | 0.11; 0.55; 0.18; 0.41 |
| 4D | displayed network | 7 metabolites, 44 genes, 62 edges | 6 metabolites, 42 genes, 60 edges: indole-3-acetic acid and its three edges (*SLC9A3*, *NOD1*, *REG3G*) are no longer displayed and propionate–*SOD2* is added; the 59 remaining edges have unchanged correlations (see `not_displayed_metabolites` below) |
| 4E | IAA–*CFB* ρ | 0.556 | 0.555 |
| S5B | module correlations | | rows of chenodeoxycholic, deoxycholic and indole-3-acetic acid change slightly |
| 5A | metabolite layer | 25 metabolites | 26 metabolites |
| 5A | KO–metabolite Procrustes | r 0.487 [0.349, 0.624], P 0.0165 | r 0.492 [0.357, 0.627], P 0.0126 |
| 5A | metabolite–host Procrustes | r 0.440 [0.315, 0.566], P 0.0925 | r 0.442 [0.315, 0.570], P 0.0840 |
| 5A | GPA agreement; mean pairwise r; P | 0.774; 0.502; 0.0007 | 0.768; 0.502; 0.0019 |
| 5B | variance explained by the metabolite view | | changes by up to 1.4 percentage points per factor (Factor 7: 3.7 % → 3.6 %) |
| 5E | Factor 7 metabolite loadings | nicotinic acid 0.350, IPA 0.251, glutamic acid 0.105 | 0.352, 0.247, 0.119 |

*MOFA.* The metabolite view changed (26 instead of 25 metabolites), so the
model grid was retrained (MOFA2 1.22.1, mofapy2 0.7.5) and the retrained
model replaces the frozen one in `input/mofa/`. The selected model again has
K = 14 factors with spike-slab weights (seed 97 instead of 41 had the highest
ELBO). Every factor of the retrained model corresponds to the factor with
the same number in the submitted model (|r| of the factor scores ≥ 0.99998);
eight factors came out with the opposite sign and were oriented like the
submitted model (`align_factor_signs()` in `R/_helpers_mofa.R`, table
`Fig5B_MOFA_sign_alignment.csv`). The statistics quoted in section 6 for
5B–F (Factor 7 P 0.011, q 0.154, max-F P 0.163; Factor 11 P 0.083; plane
R² 0.157 / P 0.020; joint R² 0.661; *TJP1* correlations) and the features of
5E and S6 are the same at the printed precision.

*Random forest (S4).* The frozen results in `input/random_forest/`
(`rf_results_submitted.rds`, `rf_inputs_frozen.rds`) still hold the
metabolite model and input of the first report (AUC 0.744), and script 09b
replays them. The metabolite model of the manuscript was refitted with the
corrected table (AUC 0.846 in Figure S4B; selected features acetate and
indole-3-lactic acid); that refit is not yet part of the frozen results.
`refit_assays <- "Metabolite"` with `freeze_refit <- TRUE` in 09b refits the
metabolite model from the corrected table and keeps the other models as
submitted.

Switches added for this re-run:

| script | switch | purpose |
|---|---|---|
| 07 | `bootstrap_sequence` | the bootstrap intervals are drawn in the metabolite order of the submitted analysis, so that metabolites whose values did not change keep their intervals; metabolites that became testable are appended |
| 08 | axis orientation after `cmdscale` | the sign of a principal coordinate is arbitrary and differed between R versions; the axes are oriented by the sign of their correlation with acetate, as in the submitted panel |
| 15 | `not_displayed_metabolites` (default `character(0)`) | histamine became quantifiable in the corrected table and enters the display candidates by the selection rule; it ends below the three-gene rule and is not drawn, and indole-3-acetic acid drops out of the panel (6 metabolites, 42 genes, 60 edges). `"histamine"` leaves it out before the ranking and gives the panel of the first report (7 metabolites, 44 genes, 62 edges) |
| 16 | `keep_submitted_row_order` | the metabolite rows keep the order that this script gave with the first report; without it the bile-acid rows are re-clustered and change places |
| 20 | `freeze_refit` | see above |

*Drawings under R 4.6.1.* With unchanged data, a few panels are not
pixel-identical to the R 4.3.2 drawings: `ggplot2` stacks several legends of
one panel in a version-dependent order (now fixed with explicit `guides()`
in 16 and 21), and label placement (`ggrepel`), point jitter and
`ComplexHeatmap` rendering differ slightly between package versions. The
underlying tables are the same.

---

## 8. Changing things

* **Colours / theme** — `R/_common.R` (`response_colors`, `family_colors`,
  `theme_tnt()`); every panel follows.
* **Thresholds** — each script keeps its thresholds in a clearly marked block
  near the top (`prevalence_cut`, `p_cut`, `n_top`, …, marked `# <--`); the
  header's WHAT CHANGES WHAT list says what moves when you change them.
* **Sample sets** — `read_stool_samples(timepoint = "Before")` etc.; the
  metabolome / host subsets are defined once in scripts 10, 14 and 18.
* **New data** — replace the files in `input/` keeping the column names in
  `input/README.md`; `dev/prepare_inputs_from_source.R` shows how the curated
  files were produced from the raw laboratory tables.

---

## 9. Citation

Please cite the manuscript when using this code or the curated inputs.
Third-party resources: MetaPhlAn 4 / HUMAnN 3 (Beghini et al. 2021), dbCAN,
KEGG Orthology, MetaCyc, MSigDB 2026.1.Hs, MOFA+ (Argelaguet et al. 2020),
vegan, mixOmics, DESeq2, fgsea, ComplexHeatmap.
