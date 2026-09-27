## =============================================================================
##  dev/prepare_inputs_from_source.R
##  ---------------------------------------------------------------------------
##  PURPOSE
##    One-time provenance script that builds the self-contained `input/`
##    directory of this repository from the authors' original data root
##    ("D:/2-연구/2-CRC metagenomics").  Everything the 28 analysis scripts
##    need is copied, de-identified and written in a plain, documented format
##    (CSV / TSV / RDS).  After this script has run once, the analysis scripts
##    never touch the original data root again.
##
##  This script is NOT part of the reproduction workflow for readers: it only
##  documents how `input/` was derived.  Readers start from `input/` directly.
##
##  DE-IDENTIFICATION RULES
##    * Trial enrolment IDs (SNU_ID) and hospital chart numbers are replaced
##      by the anonymous study SubjectID ("S22211" ...).
##    * Patient-level clinical covariates are written to a separate,
##      restricted file (input/metadata/patients_clinical.csv) that is NOT
##      meant for public upload (see input/README.md).
##
##  OUTPUT  (all under <package>/input/)
##    metadata/      patients.csv, patients_clinical.csv (restricted),
##                   stool_samples.csv, rna_samples.csv, id_crosswalk_restricted.csv
##    microbiome/    metaphlan_{sgb,species,genus,family,phylum}.csv,
##                   humann_pathway_abundance.tsv, humann_pathway_coverage.tsv,
##                   humann_ko_cpm.tsv, humann_ec_cpm_unstratified.tsv,
##                   humann_ec_2.3.1.9_species_stratified.tsv,
##                   cazy_family_tpm.tsv, cazy_subfamily_tpm.tsv,
##                   ko_functional_categories.csv, microbial_groups.csv
##    metabolome/    fecal_metabolites.csv
##    host_rnaseq/   tumor_rna_counts.csv (restricted),
##                   msigdb_2026.1.Hs_gene_sets.rds,
##                   gsea_targeted_results_frozen.csv,
##                   gsea_global_results_frozen.csv,
##                   enrichment_map_nodes_frozen.csv,
##                   fig4C_gsea_annotation.csv
##    random_forest/ rf_inputs_frozen.rds, rf_results_submitted.rds
##    mofa/          mofa_input_frozen.rds, mofa_trained_model_frozen.rds
##    INPUT_MANIFEST.csv   (file, rows, columns, md5, restricted flag)
## =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tools)
})

## ---- 0. Locate folders -------------------------------------------------------

src <- "D:/2-연구/2-CRC metagenomics"                       # original data root

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f)) return(dirname(normalizePath(f[1], winslash = "/")))
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of)) return(dirname(normalizePath(of, winslash = "/")))
  getwd()
}
pkg <- normalizePath(file.path(get_script_dir(), ".."), winslash = "/")
if (!file.exists(file.path(pkg, "LARC_TNT_response.Rproj"))) {
  pkg <- "D:/7-논문/2025 - 04 - TNT metagenomics/260922 Gut Microbes 투고본/LARC_TNT_response"
}
inp <- file.path(pkg, "input")
for (d in c("metadata", "microbiome", "metabolome", "host_rnaseq", "random_forest", "mofa")) {
  dir.create(file.path(inp, d), recursive = TRUE, showWarnings = FALSE)
}
message("Package root : ", pkg)
message("Source root  : ", src)

write_csv_utf8 <- function(x, path) {
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}
write_tsv_utf8 <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8", na = "")
}

## ---- 1. Patients and stool samples -------------------------------------------
## The expanded clinical workbook (26 baseline patients, 42 stool samples) is the
## superset of the earlier 21-column workbook; both share SampleID/TRG columns.

clin <- as.data.frame(read_xlsx(file.path(src, "input/student_final_260907/metadata_최종.xlsx")))
clin$SampleID <- paste0("Sample_", clin$SampleID)
old  <- as.data.frame(read_xlsx(file.path(src, "260224 final Input file/metadata_최종.xlsx")))
old$SampleID <- paste0("Sample_", old$SampleID)
stopifnot(setequal(clin$SampleID, old$SampleID), nrow(clin) == 42)

## Trial ID -> SubjectID crosswalk (kept only in the restricted crosswalk file)
snu_map <- clin %>% distinct(SNU_ID, SubjectID)
stopifnot(!anyDuplicated(snu_map$SNU_ID), !anyDuplicated(snu_map$SubjectID), nrow(snu_map) == 26)

## Stool sample table (public).  AnalysisOrder = row order of the original
## metadata workbook (trial ID, then visit); the original scripts processed the
## samples in this order, and permutation-based P values (PERMANOVA) depend on it.
stool_samples <- clin %>%
  transmute(AnalysisOrder = match(SampleID, old$SampleID), SampleID, SubjectID, Timepoint = TNT) %>%
  arrange(AnalysisOrder)
stopifnot(all(stool_samples$Timepoint %in% c("Before", "Ongoing")), !anyNA(stool_samples$AnalysisOrder))
message("Workbook row orders identical (clin vs old): ", identical(clin$SampleID, old$SampleID))

## EnrolmentOrder = rank of the (withheld) trial ID; an anonymous 1..26 index that
## reproduces the patient order of the original analyses
enrol <- snu_map %>% arrange(SNU_ID) %>% mutate(EnrolmentOrder = row_number()) %>% select(SubjectID, EnrolmentOrder)

## Patient table (public): response labels only
patients <- clin %>%
  filter(TNT == "Before") %>%                        # one baseline row per patient
  transmute(
    SubjectID,
    Response  = ifelse(TRG_1 == "CR", "pCR", "non-pCR"),
    TRG_score = as.integer(TRG_score),
    TRG_label = TRG
  ) %>%
  distinct() %>%
  left_join(enrol, by = "SubjectID") %>%
  arrange(EnrolmentOrder)
stopifnot(nrow(patients) == 26, !anyDuplicated(patients$SubjectID), !anyNA(patients$EnrolmentOrder))

## Patient clinical covariates (RESTRICTED)
old_pat <- old %>%
  transmute(SubjectID, cT_stage_code = as.integer(Pre_Op_Tstage)) %>%  # 2,3,4 (=cT4a),5 (=cT4b)
  distinct()
patients_clinical <- clin %>%
  filter(TNT == "Before") %>%                        # covariates recorded at the baseline visit
  transmute(
    SubjectID,
    Sex, Age = as.numeric(Age), BMI = as.numeric(BMI), ASA = as.integer(ASA),
    cT_label = cTstage, cN_label = cNstage,
    cT_stage = as.integer(Pre_Op_Tstage),           # 2/3/4  (cT4a and cT4b merged)
    cN_stage = as.integer(Pre_Op_Nstage),           # 0 / 1  (node-positive)
    AJCC_stage = as.integer(AJCCstage),
    CEA_recorded = as.character(CEA),               # as recorded ("<0.5" possible)
    CEA = as.numeric(ifelse(CEA == "<0.5", "0.5", CEA)),
    Diabetes, Hypertension, Heart_disease, Pulmonary_disease, Liver_disease,
    Cerebrovascular_accident, Smoking, Drinking, Cancer_family_history,
    Glucose = as.numeric(Glucose), Albumin = as.numeric(Albumin), Hemoglobin = as.numeric(Hemoglobin),
    CRP = as.numeric(C.reactive_protein), Segmented_neutrophil = as.numeric(Segmented_neutrophil),
    Lymphocyte = as.numeric(Lymphocyte), NLR = as.numeric(NLR), WBC = as.numeric(WBC),
    Creatinine = as.numeric(Creatine)
  ) %>%
  distinct() %>%
  left_join(old_pat, by = "SubjectID") %>%
  arrange(SubjectID)
stopifnot(nrow(patients_clinical) == 26)

write_csv_utf8(patients, file.path(inp, "metadata/patients.csv"))
write_csv_utf8(patients_clinical, file.path(inp, "metadata/patients_clinical.csv"))
write_csv_utf8(stool_samples, file.path(inp, "metadata/stool_samples.csv"))

## ---- 2. Fecal metabolome ------------------------------------------------------
## `metabolite` in metabolite_metadata.RData is the 32-sample matched table used
## by every metabolite analysis (20 baseline + 12 after-RT).

e_met <- new.env(); load(file.path(src, "input/metabolite_metadata.RData"), envir = e_met)
met <- e_met$metabolite
met_cols <- e_met$metabolite_cols
stopifnot(length(met_cols) == 40, nrow(met) == 32)

## cross-check against the raw laboratory table
raw_met <- read.csv(file.path(src, "input/metabolites.csv"), check.names = FALSE)
chk <- merge(met[, c("DNA_ID", met_cols)], raw_met[, c("SampleID", met_cols)],
             by.x = "DNA_ID", by.y = "SampleID", suffixes = c("", ".raw"))
stopifnot(nrow(chk) == 32,
          isTRUE(all.equal(as.matrix(chk[, met_cols]), as.matrix(chk[, paste0(met_cols, ".raw")]),
                           check.attributes = FALSE)))

## AnalysisOrder = row order of the laboratory table (metabolites.csv), used by
## the original ordination / PERMANOVA code
fecal_metabolites <- met %>%
  transmute(AnalysisOrder = match(DNA_ID, raw_met$SampleID), MetaboliteSampleID = DNA_ID,
            SubjectID, Timepoint = TNT, SampleID, across(all_of(met_cols))) %>%
  rename(Lithocholic_acid = Lithocholi_acid) %>%     # fix the historical typo
  arrange(AnalysisOrder)
stopifnot(!anyNA(fecal_metabolites$AnalysisOrder))
stopifnot(all(fecal_metabolites$SampleID %in% stool_samples$SampleID))
write_csv_utf8(fecal_metabolites, file.path(inp, "metabolome/fecal_metabolites.csv"))

## ---- 3. Tumor RNA-seq -----------------------------------------------------------
## RNA samples are linked to the study through hospital chart numbers; only the
## de-identified link (RNA_sample_id -> SubjectID/Timepoint) is exported.

rna_meta <- read.csv(file.path(src, "host_RNAseq/meta_RNA.csv"), check.names = FALSE)
rna_map  <- read.csv(file.path(src, "host_RNAseq/chart_numb_matching.csv"), check.names = FALSE)
rna_link <- clin %>%
  transmute(SampleID, SubjectID, Timepoint = TNT, SNU_AL_ID = SNU_ID) %>%
  inner_join(rna_map %>% select(Chart_numb, SNU_AL_ID), by = "SNU_AL_ID") %>%
  mutate(key = paste0(Chart_numb, "_", ifelse(Timepoint == "Before", 1, 2))) %>%
  inner_join(
    rna_meta %>%
      filter(tolower(timepoint) %in% c("pre", "post")) %>%
      transmute(RNA_sample_id, key = paste0(Chart_numb, "_", ifelse(tolower(timepoint) == "post", 2, 1))),
    by = "key"
  ) %>%
  distinct(RNA_sample_id, .keep_all = TRUE)

rna_samples <- rna_link %>%
  transmute(RNA_sample_id, SubjectID, Timepoint, Tissue = "tumor", SampleID) %>%
  arrange(SubjectID, Timepoint)
stopifnot(nrow(rna_samples) == 27, !anyDuplicated(rna_samples$RNA_sample_id),
          sum(rna_samples$Timepoint == "Before") == 20)
write_csv_utf8(rna_samples, file.path(inp, "metadata/rna_samples.csv"))

rna <- read.csv(file.path(src, "host_RNAseq/TNT_Expression_Profile.GRCh38.gene.csv"), check.names = FALSE)
count_cols <- paste0(rna_samples$RNA_sample_id, "_Read_Count")
stopifnot(all(count_cols %in% names(rna)))
tumor_counts <- rna[, c("Gene_ID", "Gene_Symbol", "Transcript_ID", "gene_biotype", count_cols)]
names(tumor_counts)[-(1:4)] <- rna_samples$RNA_sample_id
write_csv_utf8(tumor_counts, file.path(inp, "host_rnaseq/tumor_rna_counts.csv"))

## Restricted crosswalk (never upload): lets the authors trace IDs back
write_csv_utf8(
  rna_link %>% transmute(RNA_sample_id, SampleID, SubjectID, SNU_ID = SNU_AL_ID, Chart_numb),
  file.path(inp, "metadata/id_crosswalk_restricted.csv")
)

## ---- 4. MetaPhlAn taxonomic profiles -----------------------------------------

for (rank in c("strain", "species", "genus", "family", "phylum")) {
  x <- read.csv(file.path(src, "260224 final Input file/Abundance file", paste0(rank, ".csv")), check.names = FALSE)
  stopifnot(setequal(setdiff(names(x), "clade_name"), stool_samples$SampleID))
  x <- x[, c("clade_name", stool_samples$SampleID)]
  out <- if (rank == "strain") "metaphlan_sgb.csv" else paste0("metaphlan_", rank, ".csv")
  write_csv_utf8(x, file.path(inp, "microbiome", out))
}

## ---- 5. HUMAnN pathways, KOs and ECs -------------------------------------------

pa <- read.delim(file.path(src, "260224 final Input file/Abundance file/merged_pathabundance.tsv"), check.names = FALSE, quote = "")
pc <- read.delim(file.path(src, "260224 final Input file/Abundance file/merged_pathcoverage.tsv"), check.names = FALSE, quote = "")
names(pa)[1] <- names(pc)[1] <- "Pathway"
stopifnot(!anyDuplicated(pa$Pathway), setequal(pa$Pathway, pc$Pathway), setequal(names(pa)[-1], stool_samples$SampleID))
pc <- pc[match(pa$Pathway, pc$Pathway), ]          # same row order as the abundance table
write_tsv_utf8(pa[, c("Pathway", stool_samples$SampleID)], file.path(inp, "microbiome/humann_pathway_abundance.tsv"))
write_tsv_utf8(pc[, c("Pathway", stool_samples$SampleID)], file.path(inp, "microbiome/humann_pathway_coverage.tsv"))

ko <- read.delim(file.path(src, "Script (260824)/Input/merged_genefamilies_KO_named.tsv"), check.names = FALSE, quote = "")
names(ko)[1] <- "GeneFamily"
ko <- ko[!grepl("|", ko$GeneFamily, fixed = TRUE), ]                       # unstratified rows only
ko <- ko[grepl("^K[0-9]{5}", ko$GeneFamily), ]
ko_out <- data.frame(
  KO      = sub(":.*$", "", ko$GeneFamily),
  KO_name = sub("^K[0-9]{5}:\\s*", "", ko$GeneFamily),
  ko[, stool_samples$SampleID], check.names = FALSE
)
stopifnot(!anyDuplicated(ko_out$KO))
write_tsv_utf8(ko_out, file.path(inp, "microbiome/humann_ko_cpm.tsv"))

ec <- read.delim(file.path(src, "Script (260824)/Input/merged_genefamilies_EC_named.tsv"), check.names = FALSE, quote = "")
names(ec)[1] <- "GeneFamily"
ec_unstrat <- ec[!grepl("|", ec$GeneFamily, fixed = TRUE), ]
ec_unstrat <- data.frame(
  EC      = sub(":.*$", "", ec_unstrat$GeneFamily),
  EC_name = sub("^[^:]*:\\s*", "", ec_unstrat$GeneFamily),
  ec_unstrat[, stool_samples$SampleID], check.names = FALSE
)
write_tsv_utf8(ec_unstrat, file.path(inp, "microbiome/humann_ec_cpm_unstratified.tsv"))
ec_2319 <- ec[grepl("^2\\.3\\.1\\.9:", ec$GeneFamily) & grepl("|", ec$GeneFamily, fixed = TRUE), ]
ec_2319 <- data.frame(
  EC = "2.3.1.9", Taxon = sub("^[^|]*\\|", "", ec_2319$GeneFamily),
  ec_2319[, stool_samples$SampleID], check.names = FALSE
)
write_tsv_utf8(ec_2319, file.path(inp, "microbiome/humann_ec_2.3.1.9_species_stratified.tsv"))

## ---- 6. CAZymes ------------------------------------------------------------------

caz_fam <- read.delim(file.path(src, "260224 final Input file/CAZyme/CAZy_family_TPM_matrix.tsv"), check.names = FALSE)
caz_sub <- read.delim(file.path(src, "260224 final Input file/CAZyme/CAZy_subfamily_TPM_matrix.tsv"), check.names = FALSE)
names(caz_fam)[1] <- "CAZy_family"; names(caz_sub)[1] <- "CAZy_subfamily"
write_tsv_utf8(caz_fam[, c("CAZy_family", stool_samples$SampleID)], file.path(inp, "microbiome/cazy_family_tpm.tsv"))
write_tsv_utf8(caz_sub[, c("CAZy_subfamily", stool_samples$SampleID)], file.path(inp, "microbiome/cazy_subfamily_tpm.tsv"))

## ---- 7. KO functional categories (frozen annotation; Supplementary Table 3) ----

kc <- readRDS(file.path(src, "input/260619 KO_category_v3.RDS"))
ko_categories <- data.frame(
  KO = kc$KO, KO_name = kc$name,
  upper_category = as.character(kc$upper_category),
  assigned_category = as.character(kc$assigned_category),
  upper_category_order = as.integer(kc$upper_category),
  assigned_category_order = as.integer(kc$assigned_category),
  confidence = kc$confidence, brite_category = kc$brite_category,
  stringsAsFactors = FALSE
)
stopifnot(nrow(ko_categories) == 2228)
write_csv_utf8(ko_categories, file.path(inp, "microbiome/ko_functional_categories.csv"))

## ---- 8. Curated microbial groups (Supplementary Table 1, as used by the code) ----

grp <- list(
  health_S1A = c("GGB13404_SGB14252", "GGB4552_SGB6276", "Clostridium_sp_NSJ_42", "Faecalibacterium_prausnitzii",
                 "GGB4585_SGB6340", "GGB9237_SGB14179", "GGB9707_SGB15229", "Intestinimonas_gabonensis",
                 "GGB3653_SGB4964", "GGB3619_SGB4894", "Lachnospiraceae_bacterium_OM04_12BH",
                 "Oscillibacter_sp_MSJ_31", "Roseburia_sp_BX1005"),
  crc_S1A = c("Dialister_pneumosintes", "Fusobacterium_nucleatum", "Fusobacterium_periodonticum",
              "Gemella_morbillorum", "Granulicatella_adiacens", "Parvimonas_micra",
              "Peptostreptococcus_stomatis", "Prevotella_intermedia", "Prevotella_nigrescens",
              "Slackia_exigua", "Solobacterium_moorei", "Solobacterium_SGB6833", "Streptococcus_mutans"),
  icb_S1A = c("Roseburia_sp_AF02_12", "Dorea_formicigenerans", "Wujia_chipingensis", "Phocaeicola_dorei",
              "Anaerofilum_hominis", "Lacrimispora_amygdalina", "Clostridium_sp_AM22_11AC", "GGB9699_SGB15216",
              "Alistipes_communis", "Anaeromassilibacillus_senegalensis", "GGB47687_SGB2286", "GGB9634_SGB15093",
              "Holdemania_filiformis", "Bacteroides_stercoris", "Blautia_caecimuris"),
  mucin_S1A = c("Akkermansia_muciniphila", "Bacteroides_caccae", "Bacteroides_clarus", "Bacteroides_faecis",
                "Bacteroides_finegoldii", "Bacteroides_fragilis", "Bacteroides_intestinalis", "Bacteroides_nordii",
                "Bacteroides_ovatus", "Bacteroides_stercoris", "Bacteroides_thetaiotaomicron", "Bacteroides_uniformis",
                "Bacteroides_xylanisolvens", "Phocaeicola_vulgatus", "Barnesiella_intestinihominis",
                "Ruminococcus_gnavus", "Ruminococcus_torques", "Bifidobacterium_bifidum", "Bifidobacterium_breve",
                "Bifidobacterium_longum"),
  butyrate_S1A = c("Anaerobutyricum_hallii", "Anaerobutyricum_soehngenii", "Anaerostipes_butyraticus",
                   "Anaerostipes_caccae", "Anaerostipes_hadrus", "Anaerotruncus_colihominis", "Coprococcus_comes",
                   "Coprococcus_eutactus", "Eubacterium_limosum", "Eubacterium_ramulus", "Eubacterium_rectale",
                   "Faecalibacterium_prausnitzii", "Roseburia_faecis", "Roseburia_hominis", "Roseburia_intestinalis",
                   "Ruminococcus_bromii", "Ruminococcus_callidus", "Coprococcus_catus", "Roseburia_inulinivorans"),
  sulfidogenic_S1A = c("Desulfovibrio_desulfuricans", "Desulfovibrio_piger", "Desulfovibrio_porci", "Bilophila_wadsworthia"),
  butyrate_producer_S1B = c("Eubacterium_rectale", "Eisenbergiella_tayi", "Roseburia_inulinivorans", "Intestinimonas_butyriciproducens"),
  butyrate_crossfeeder_S1B = c("Bifidobacterium_longum", "Ruminococcus_callidus", "Bifidobacterium_dentium", "Lacrimispora_amygdalina"),
  icb_S1B = c("Eubacterium_rectale", "Bifidobacterium_longum", "Ruminococcus_callidus", "Bifidobacterium_dentium",
              "Lachnospiraceae_bacterium_CLA_AA_H244", "Akkermansia_muciniphila", "Lacrimispora_amygdalina",
              "Ruminococcus_sp_AF41_9"),
  crc_S1B = c("Veillonella_dispar", "Akkermansia_muciniphila", "Schaalia_odontolytica", "Clostridium_saudiense",
              "Turicibacter_sanguinis", "Rothia_mucilaginosa", "Schaalia_SGB17154", "Solobacterium_moorei",
              "Lancefieldella_parvula", "Fusobacterium_nucleatum")
)
group_label <- c(
  health_S1A = "Cardiometabolic health-associated bacteria",
  crc_S1A = "Oral-derived CRC signatures", icb_S1A = "ICB-favorable bacteria",
  mucin_S1A = "Mucin degraders", butyrate_S1A = "Butyrate-producing guild",
  sulfidogenic_S1A = "Sulfidogenic bacteria",
  butyrate_producer_S1B = "Direct butyrate producer", butyrate_crossfeeder_S1B = "Butyrate cross-feeder",
  icb_S1B = "ICB-favorable bacteria", crc_S1B = "CRC signatures"
)
microbial_groups <- bind_rows(lapply(names(grp), function(k) {
  data.frame(List = k, Table = ifelse(grepl("S1A$", k), "Supplementary Table 1A", "Supplementary Table 1B"),
             Group = group_label[[k]], Species = unique(grp[[k]]), stringsAsFactors = FALSE)
}))
sp_names <- sub(".*s__", "", read.csv(file.path(inp, "microbiome/metaphlan_species.csv"), check.names = FALSE)$clade_name)
stopifnot(all(microbial_groups$Species %in% sp_names))
write_csv_utf8(microbial_groups, file.path(inp, "microbiome/microbial_groups.csv"))

## ---- 9. Frozen host enrichment results ------------------------------------------

e_gs <- new.env()
load(file.path(src, "host_RNAseq/results_clean_metabolite_host/pathway_metabolite/host_metabolite_pathway_analysis_results.RData"), envir = e_gs)
saveRDS(
  list(
    gene_sets = e_gs$gene_sets_use,                     # named list of 2,736 gene sets (MSigDB 2026.1.Hs)
    gene_set_table = e_gs$gene_set_size_table,          # collection / description / size
    gene_sets_long = e_gs$gene_sets_use_long,           # long gene-membership table
    msigdb_version = "2026.1.Hs"
  ),
  file.path(inp, "host_rnaseq/msigdb_2026.1.Hs_gene_sets.rds"), compress = "xz"
)
fg <- e_gs$fgsea_res
fg$leadingEdge <- vapply(fg$leadingEdge, function(z) paste(unlist(z), collapse = ";"), character(1))
write_csv_utf8(fg, file.path(inp, "host_rnaseq/gsea_targeted_results_frozen.csv"))

arch <- file.path(src, "host_RNAseq/results_clean_metabolite_host/Figure5_host_RNAseq")
file.copy(file.path(arch, "Fig5A_unified_GSEA_results_with_direction.csv"),
          file.path(inp, "host_rnaseq/gsea_global_results_frozen.csv"), overwrite = TRUE)
file.copy(file.path(arch, "Fig5A_enrichment_map_nodes_mixed_direction.csv"),
          file.path(inp, "host_rnaseq/enrichment_map_nodes_frozen.csv"), overwrite = TRUE)

e_f5 <- new.env(); load(file.path(arch, "Fig5_host_RNAseq_inputs.RData"), envir = e_f5)
gpd <- as.data.frame(e_f5$gsea_plot_df)
gpd$leadingEdge <- as.character(gpd$leadingEdge)
write_csv_utf8(
  gpd[, c("main_category", "family_id", "display_label", "Host_axis", "pathway_label", "pathway_id",
          "gs_name", "NES", "pval", "padj", "size", "leadingEdge", "pass_display_filter")],
  file.path(inp, "host_rnaseq/fig4C_gsea_annotation.csv")
)

## ---- 10. Random forest: frozen inputs and submitted results ---------------------

e_rf <- new.env(parent = emptyenv())
load(file.path(src, "input/student_final_260907/05_RF modeling 260907.RData"), envir = e_rf)
rf_meta <- e_rf$m
rf_meta$SubjectID <- snu_map$SubjectID[match(rf_meta$SNU_ID, snu_map$SNU_ID)]
stopifnot(!anyNA(rf_meta$SubjectID))
assay_objects <- c(Strain = "strain", Species = "species", Genus = "genus", Family = "family",
                   Phylum = "phylum", Pathway = "path", KO = "ko", Metabolite = "met")
rf_inputs <- lapply(names(assay_objects), function(nm) {
  prefix <- assay_objects[[nm]]
  values <- get(paste0(prefix, "_mat"), envir = e_rf)
  saved  <- get(paste0("RF_", prefix), envir = e_rf)
  sm <- rf_meta[match(rownames(values), rf_meta$SampleID), ]
  stopifnot(identical(rownames(values), sm$SampleID), identical(rownames(values), rownames(saved)),
            isTRUE(all.equal(unname(values), unname(as.matrix(saved[, colnames(values)])))),
            identical(as.character(sm$SNU_ID), as.character(saved$SNU_ID)))
  list(values = values,
       metadata = data.frame(SampleID = sm$SampleID, SubjectID = sm$SubjectID, Timepoint = sm$TNT,
                             Response = ifelse(sm$TRG_1 == "CR", "pCR", "non-pCR"), stringsAsFactors = FALSE))
})
names(rf_inputs) <- names(assay_objects)
rf_clinical <- e_rf$RF_ci
rf_clinical <- data.frame(
  SubjectID = snu_map$SubjectID[match(rf_clinical$SNU_ID, snu_map$SNU_ID)],
  Response  = ifelse(as.character(rf_clinical$TRG_1) == "CR", "pCR", "non-pCR"),
  cT_stage_code = as.integer(rf_clinical$Pre_Op_Tstage),
  cN_stage = as.integer(rf_clinical$Pre_Op_Nstage),
  CEA_high = as.integer(rf_clinical$CEA_bin), stringsAsFactors = FALSE
)
stopifnot(nrow(rf_clinical) == 26, !anyNA(rf_clinical$SubjectID))
saveRDS(list(inputs = rf_inputs, clinical = rf_clinical,
             note = "Exact assay matrices used for the submitted random-forest models (student workspace 2026-09-07)."),
        file.path(inp, "random_forest/rf_inputs_frozen.rds"), compress = "xz")

rf_res <- readRDS(file.path(src, "input/student_final_260907/260907 RF_Global_Direction_all_results.rds"))
for (nm in names(rf_res)) {
  p <- rf_res[[nm]]$Pred
  p$SubjectID <- snu_map$SubjectID[match(as.character(p$SNU_ID), snu_map$SNU_ID)]
  stopifnot(!anyNA(p$SubjectID))
  p$SNU_ID <- NULL
  if (nm == "CI") rownames(p) <- p$SubjectID
  rf_res[[nm]]$Pred <- p
}
saveRDS(rf_res, file.path(inp, "random_forest/rf_results_submitted.rds"), compress = "xz")

## ---- 11. MOFA: frozen input and trained-model tables ------------------------------

mofa_in <- readRDS(file.path(src, "figures/mofa/v66.1/results/condition_input.rds"))
keep_meta <- c("SampleID", "SubjectID", "Timepoint", "TRG_plot", "TRG_score", "Response", "Time", "n_observed_views")
mofa_in$metadata <- mofa_in$metadata[, keep_meta]
saveRDS(mofa_in, file.path(inp, "mofa/mofa_input_frozen.rds"), compress = "xz")

mofa_tr <- readRDS(file.path(src, "figures/mofa/v66.1/results/trained_model_tables.rds"))
mofa_tr$extracted$scores <- mofa_tr$extracted$scores[, c(keep_meta, grep("^Factor", names(mofa_tr$extracted$scores), value = TRUE))]
mofa_tr$training_log$path <- basename(mofa_tr$training_log$path)
mofa_tr$selected_hdf5 <- basename(mofa_tr$selected_hdf5)
saveRDS(mofa_tr, file.path(inp, "mofa/mofa_trained_model_frozen.rds"), compress = "xz")

## ---- 12. Legacy 10-patient species subset used by Figure 4E ---------------------

e_co <- new.env(); load(file.path(src, "input/coherence_data.RData"), envir = e_co)
paired_subjects <- e_co$coherence_data$paired_subjects
met_pairs <- fecal_metabolites %>% count(SubjectID) %>% filter(n == 2) %>% pull(SubjectID)
message("Legacy paired subjects equal metabolome-paired subjects: ", setequal(paired_subjects, met_pairs))
write_csv_utf8(data.frame(SubjectID = sort(paired_subjects)), file.path(inp, "metadata/legacy_paired_subjects.csv"))

## ---- 13. Manifest ------------------------------------------------------------------

restricted <- c("metadata/patients_clinical.csv", "metadata/id_crosswalk_restricted.csv",
                "host_rnaseq/tumor_rna_counts.csv")
files <- list.files(inp, recursive = TRUE)
files <- files[!files %in% c("INPUT_MANIFEST.csv", "README.md")]
manifest <- bind_rows(lapply(files, function(f) {
  p <- file.path(inp, f)
  dims <- tryCatch({
    if (grepl("\\.csv$", f)) { x <- read.csv(p, check.names = FALSE); c(nrow(x), ncol(x)) }
    else if (grepl("\\.tsv$", f)) { x <- read.delim(p, check.names = FALSE, nrows = 1); c(length(readLines(p)) - 1, ncol(x)) }
    else c(NA, NA)
  }, error = function(e) c(NA, NA))
  data.frame(file = f, rows = dims[1], columns = dims[2], size_MB = round(file.size(p) / 1e6, 2),
             md5 = unname(md5sum(p)), restricted = f %in% restricted, stringsAsFactors = FALSE)
}))
write_csv_utf8(manifest, file.path(inp, "INPUT_MANIFEST.csv"))
print(manifest)
message("DONE: input directory written to ", inp)
