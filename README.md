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
| `19_Fig5A_Procrustes_GPA.R` | 5A | pairwise Procrustes / protest and generalized Procrustes analysis with a patient-level permutation test | cache of 18 |
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
you want to reproduce permutation P values to the last digit.

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
* **metabolome/** — 40 targeted fecal metabolites (µmol/g) for 32 samples.
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
| MOFA (5B–F, S6) | `mofa/mofa_input_frozen.rds`, `mofa_trained_model_frozen.rds` | the four views are rebuilt from the raw data and must equal the frozen training input (tolerance 1e-7); variance explained is reconstructed from scores × weights and compared with the stored values | `refit_mofa <- TRUE` in 20 (needs MOFA2 + Python) |

All other analyses are recomputed from the raw tables every time.

---

## 6. Verification against the submitted figures

All 28 scripts were run from a clean `output/` folder (R 4.3.2, Windows).
Every numeric result printed in a panel or quoted in the text was checked:

* diversity P values (evenness 0.047, inverse Simpson 0.087, Shannon 0.148,
  richness 0.775), PERMANOVA S1D (R² 0.044, P 0.291), S1E R² values;
* PLS-DA VIPs (1H), KO category Z / P (2B), KO stars (2C), CAZyme P values
  (2D–F), MetaCyc pathway P values and group tests (S2);
* Hedges' g and bootstrap intervals of all seven metabolites of 3B, ratio
  test P = 0.023 (3F), metabolome PERMANOVA (R² 0.134, P 0.063), *E. rectale*
  EC 2.3.1.9 rho 0.52 / P 0.020 (3D), bile-acid index P values (S3B);
* random-forest AUCs (S4B), the 22 gene-set statistics (4B), Fig 4C stars,
  the 4D / 4E node and edge sets;
* pairwise Procrustes r 0.715 / 0.487 / 0.440 with their intervals and P
  values, GPA agreement 0.774, mean pairwise r 0.502, P 0.0007 (5A); MOFA
  Factor 7 P 0.011 (BH q 0.154, max-F P 0.163), Factor 11 P 0.083, plane
  R² 0.157 / P 0.020 (5B–D); joint R² 0.661; the 12 Factor 7 features and
  their loadings (5E); TJP1 Pearson r 0.48, Spearman 0.60, partial 0.62 (5F).

Permutation P values are reproduced to the last printed digit wherever the
sample order of the original run is known (see `AnalysisOrder` /
`EnrolmentOrder`, section 4). The three Figure 5F permutation P values
(999 permutations) agree with the submission within ±0.002.

Known, documented differences from the submitted PDF panels:

| panel | difference | where explained |
|---|---|---|
| 1A, 3E, S4A | the timeline, the tryptophan-pathway scheme and the workflow scheme are drawn with plain boxes; the submitted versions were finished in Illustrator | header of 00 |
| S1C | the submitted legend lists three families (Enterobacteriaceae, Lactobacillaceae, Selenomonadaceae) that are not in the baseline top-10 (legend artefact of the after-RT panel); reproduced optionally | 02, `s1c_legend_extra` |
| 2C | two KOs (K13542, K20491) show `**` here where the submitted panel shows `***` — P just around 0.001 depending on the Wilcoxon variant | header of 05 |
| S2A | the colour bar is labelled log10 ratio (the submitted legend said Log2FC; the values are log10) | header of 06b |
| 1G, 4A, S6 | network node positions come from a seeded force-directed layout; the submitted panels were arranged by hand (1G reproduces the Cytoscape coordinates where available) | 03, 12, 22 |
| S5B | lithocholic acid is excluded explicitly: in the laboratory table used for the submission its column was misspelled (`Lithocholi_acid`) and therefore not recognised as a bile acid; the curated input corrects the name | 16, `excluded_metabolites` |
| 4B | the first gene set ("PI3K–AKT signaling with NF-κB") is the KEGG MEDICUS set N01339 (NNK/NNN → PI3K signaling); the display name is the submitted one | 11, `catalog` |
| 5B–F, S6 | the rebuilt MOFA views equal the frozen training input (species, metabolite, host exactly; KO within 1.5 × 10⁻⁵ relative, the same 2,000 features); the frozen model is used, so a retrain (`refit_mofa`) will give different factor numbers | 20 |
| 5E | the lollipop shows the normalised loading (sign oriented towards pCR) without the Hedges' g colour of the original draft script; g is in the table | 22 |

Design: colours, shapes, axis labels and panel proportions follow the
submitted PDFs (colour codes were read from the PDFs). Fonts differ slightly
(the PDFs use Arial; R uses the system sans-serif) and a few labels in dense
panels are placed by `ggrepel` and may land a few millimetres from where they
were in the submission.

---

## 7. Changing things

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

## 8. Citation

Please cite the manuscript when using this code or the curated inputs.
Third-party resources: MetaPhlAn 4 / HUMAnN 3 (Beghini et al. 2021), dbCAN,
KEGG Orthology, MetaCyc, MSigDB 2026.1.Hs, MOFA+ (Argelaguet et al. 2020),
vegan, mixOmics, DESeq2, fgsea, ComplexHeatmap.
