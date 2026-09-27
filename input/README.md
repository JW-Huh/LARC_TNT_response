# input/ — curated input files

All files were exported once from the laboratory tables by
`dev/prepare_inputs_from_source.R`; `INPUT_MANIFEST.csv` lists every file with
its row / column count, size and MD5 checksum so that a downloaded copy can be
verified. Files marked **restricted** contain patient-level clinical or identifying
data and are not part of the public repository (see the main README for what
still runs without them).

Conventions

* CSV files are UTF-8 with a header row; TSV files are tab-separated without
  quotes. Missing values are empty cells (read as `NA` by
  `read_csv_na()` in `R/_common.R`).
* Sample columns of the microbiome tables are named `Sample_1` … `Sample_42`
  (= `SampleID` in `metadata/stool_samples.csv`).
* Identifiers are anonymised: `SubjectID` (`S22207` …), stool `SampleID`
  (`Sample_7` …), `RNA_sample_id` (`TNT_RNA_44` …), `MetaboliteSampleID`
  (laboratory metabolome sample code).
* Taxon names use the MetaPhlAn 4 clade spelling (`Eubacterium_rectale`,
  `Ruminococcus_sp_AF13_28`); `pretty_species()` / `taxon_plotmath()` in
  `R/_common.R` convert them for display.

---

## metadata/

| file | rows | columns | notes |
|---|---|---|---|
| `patients.csv` | 26 | `SubjectID`, `Response` (`pCR` / `non-pCR`), `TRG_score` (0–3), `TRG_label` (`CR` = 0, `nearCR` = 1, `PR` = 2, `Poor` = 3), `EnrolmentOrder` (1–26) | one row per patient. `Response = pCR` ⇔ `TRG_score = 0`. `EnrolmentOrder` is an anonymous rank that reproduces the patient order of the original analysis (bootstrap / permutation results depend on it) |
| `stool_samples.csv` | 42 | `AnalysisOrder`, `SampleID`, `SubjectID`, `Timepoint` (`Before` = baseline, `Ongoing` = after radiotherapy) | 26 baseline + 16 after-RT stool metagenomes. `AnalysisOrder` = the sample order of the original metadata workbook (keep it: PERMANOVA / envfit P values depend on row order) |
| `rna_samples.csv` | 27 | `RNA_sample_id`, `SubjectID`, `Timepoint`, `Tissue` (`tumor`), `SampleID` (matching stool sample) | 20 baseline + 7 after-RT tumour biopsies |
| `legacy_paired_subjects.csv` | 10 | `SubjectID` | the ten patients with both a baseline and an after-RT metabolome; the species edges of Figure 4E were computed on this subset in the original analysis (kept for exact reproduction) |
| `patients_covariates.csv` | 26 | `SubjectID`, `Sex` (`M`/`F`), `Age` (years), `BMI` | the three covariates used by the adjusted correlations (S5B, S6, Figure 5F). Public, de-identified minimum of the clinical table |
| `sfig1E_permanova_frozen.csv` | 26 | `Covariate`, `N`, `Type`, `R2_percent`, `F`, `P`, `Group` | per-covariate PERMANOVA results of the submitted run (Supplementary Figure 1E); lets script 04b draw the panel without the restricted clinical table |
| `patients_clinical.csv` **restricted** | 26 | `SubjectID`, `Sex` (`M`/`F`), `Age`, `BMI`, `ASA`, `cT_label`, `cN_label`, `cT_stage`, `cN_stage`, `AJCC_stage`, `CEA_recorded`, `CEA`, `Diabetes`, `Hypertension`, `Heart_disease`, `Pulmonary_disease`, `Liver_disease`, `Cerebrovascular_accident`, `Smoking`, `Drinking`, `Cancer_family_history`, `Glucose`, `Albumin`, `Hemoglobin`, `CRP`, `Segmented_neutrophil`, `Lymphocyte`, `NLR`, `WBC`, `Creatinine`, `cT_stage_code` | binary history variables are coded `O` / `X`; smoking / drinking `Present` / `Past` / `X`; empty cells = not recorded. Used by S1E (all columns), S5B, S6 and 5F (Age, Sex, BMI only) |
| `id_crosswalk_restricted.csv` **restricted** | 27 | `RNA_sample_id`, `SampleID`, `SubjectID`, `SNU_ID`, `Chart_numb` | laboratory identifiers; only needed to trace samples back to the raw data |

## microbiome/

| file | rows × sample cols | id columns | content |
|---|---|---|---|
| `metaphlan_sgb.csv` | 1,317 × 42 | `clade_name` (full MetaPhlAn lineage ending in `t__SGBxxxx`) | MetaPhlAn 4 (vJun23 database) relative abundance (%) at species-level genome bin (SGB) resolution; columns sum to 100 |
| `metaphlan_species.csv` | 1,227 × 42 | `clade_name` (… `s__Genus_species`) | species-level profiles |
| `metaphlan_genus.csv`, `metaphlan_family.csv`, `metaphlan_phylum.csv` | 657 / 241 / 13 × 42 | `clade_name` | genus / family / phylum profiles |
| `humann_ko_cpm.tsv` | 7,205 × 42 | `KO`, `KO_name` | HUMAnN 3 gene families regrouped to KEGG Orthology, unstratified, copies per million (CPM); rows `UNMAPPED`, `UNGROUPED` and ribosomal proteins are removed by the scripts where stated |
| `humann_ec_cpm_unstratified.tsv` | 2,554 × 42 | `EC`, `EC_name` | HUMAnN 3 gene families regrouped to Enzyme Commission numbers (CPM) |
| `humann_ec_2.3.1.9_species_stratified.tsv` | 106 × 42 | `EC` (2.3.1.9), `Taxon` (`g__…\|s__…`) | contribution of each species to acetyl-CoA C-acetyltransferase (CPM); Figure 3D |
| `humann_pathway_abundance.tsv` | 23,748 × 42 | `Pathway` (`PWY-621: sucrose degradation III…` or `PWY-621\|g__Bacteroides.s__…` for stratified rows) | HUMAnN 3 MetaCyc pathway abundance (unstratified + taxon-stratified rows) |
| `humann_pathway_coverage.tsv` | 23,748 × 42 | `Pathway` | matching MetaCyc pathway coverage (0–1); a pathway is masked to 0 where coverage < 0.5 |
| `cazy_family_tpm.tsv` | 366 × 42 | `CAZy_family` (`GH13`, `CBM20` …) | dbCAN annotation of assembled genes, transcripts per million |
| `cazy_subfamily_tpm.tsv` | 6,528 × 42 | `CAZy_subfamily` (`GH13_31` …) | subfamily-level TPM |
| `ko_functional_categories.csv` | 2,228 | `KO`, `KO_name`, `upper_category` (5), `assigned_category` (19), `upper_category_order`, `assigned_category_order`, `confidence`, `brite_category` | rule-based KO → functional category assignment (Supplementary Table 3); 2,127 unique KOs, 101 KOs appear in two categories |
| `microbial_groups.csv` | 110 | `List`, `Table`, `Group`, `Species` | curated group memberships of Supplementary Table 1: `health_S1A`, `crc_S1A`, `crc_S1B`, `icb_S1A`, `icb_S1B`, `mucin_S1A`, `butyrate_S1A`, `butyrate_producer_S1B`, `butyrate_crossfeeder_S1B`, `sulfidogenic_S1A` (Figure 1E–F annotations) |

## metabolome/

| file | rows | columns | content |
|---|---|---|---|
| `fecal_metabolites.csv` | 32 | `AnalysisOrder`, `MetaboliteSampleID`, `SubjectID`, `Timepoint`, `SampleID`, then 40 metabolite columns | targeted fecal metabolomics (µmol per g stool): 6 short-chain fatty acids, 18 tryptophan-pathway / amine / amino-acid analytes, 16 bile acids. 20 baseline + 12 after-RT samples (21 patients). Values at the assay floor are repeated identical minima (see the QC steps in scripts 07, 14, 18). `AnalysisOrder` = the row order of the laboratory table (keep it). `Lithocholic_acid` was spelled `Lithocholi_acid` in the laboratory table (see script 16) |

## host_rnaseq/

| file | content |
|---|---|
| `tumor_rna_counts.csv` | gene-level read counts of the 27 tumour biopsies (the matrix deposited in GEO GSE339362): `Gene_ID`, `Gene_Symbol`, `Transcript_ID`, `gene_biotype`, one column per `RNA_sample_id` (46,425 genes) |
| `gsea_targeted_results_frozen.csv` | fgsea results (2,803 rows = gene set × host axis) of the targeted GSEA behind Figure 4B: `pathway`, `pval`, `padj`, `log2err`, `ES`, `NES`, `size`, `leadingEdge` (`;`-separated genes), plus the gene-set annotation columns |
| `gsea_global_results_frozen.csv` | fgsea results of the whole gene-set family (8,240 sets) behind Figure 4A, with `Direction`, `display_label`, `biological_axis`, `leadingEdge_string` |
| `enrichment_map_nodes_frozen.csv` | the 160 nodes of the submitted enrichment map with their Louvain cluster (`network_cluster`), cluster label and direction |
| `fig4C_gsea_annotation.csv` | 25 gene sets with the family / category annotation used by Figure 4C (`main_category`, `family_id`, `display_label`, `gs_name`, `NES`, `pval`, `padj`, `leadingEdge`) |
| `msigdb_2026.1.Hs_gene_sets.rds` | list with `gene_sets` (named list of 2,736 MSigDB 2026.1.Hs gene sets), `gene_set_table`, `gene_sets_long`; needed only for a GSEA refit |

## random_forest/

| file | content |
|---|---|
| `rf_inputs_frozen.rds` | `inputs[[assay]]$values` (samples × features) and `$metadata` (`SampleID`, `SubjectID`, `Timepoint`, `Response`) for the eight assays (SGB, species, genus, family, phylum, MetaCyc pathway, KO, metabolite) exactly as used to train the submitted models; `clinical` (26 patients: `cT_stage_code`, `cN_stage`, `CEA_high`) |
| `rf_results_submitted.rds` | per model: selected features, out-of-fold predictions (`Pred`), AUC, importance (MDA, MDG), SHAP summaries and stability; replayed by scripts 09b / 09c |

## mofa/

| file | content |
|---|---|
| `mofa_input_frozen.rds` | `values` (four centred views: species 392, KO 2,000, metabolite 25, host 2,000 features; 42 samples, NA where a view is missing) and `metadata` (`SampleID`, `SubjectID`, `Timepoint`, `TRG_plot`, `TRG_score`, `Response`, `Time`, `n_observed_views`) — the exact training input of the frozen model; script 20 rebuilds it from the raw tables and checks equality |
| `mofa_trained_model_frozen.rds` | `extracted$scores` (sample × 14 factors), `extracted$weights` (view, feature, factor, weight, feature_label), `extracted$variance` (r² per view × factor), `training_log` (K, spike-slab, seed, ELBO, selected) of the 4 × 2 model grid |

---

### Re-creating the inputs from the raw data

`dev/prepare_inputs_from_source.R` reads the laboratory tables (path in
`src`), anonymises the identifiers, writes every file above and the manifest.
It is kept for provenance; running it requires the raw data and the
restricted crosswalk.
