## =============================================================================
##  09a_SFig4_Random_forest_inputs.R
##  ---------------------------------------------------------------------------
##  PANELS      (no figure) - prepares the input matrices for scripts 09b / 09c
##              (Supplementary Figure 4B-J)
##
##  PURPOSE     Load and check the eight assay matrices (five taxonomic ranks,
##              MetaCyc pathways, KOs, metabolites) and the clinical table that
##              were used to train the submitted random-forest models, and
##              store them in output/cache/ for the modelling script.
##
##  INPUT       input/random_forest/rf_inputs_frozen.rds
##                $inputs[[assay]]$values     samples x features matrix
##                $inputs[[assay]]$metadata   SampleID, SubjectID, Timepoint, Response
##                $clinical                   SubjectID, Response, cT_stage_code,
##                                            cN_stage, CEA_high  (26 patients)
##
##  OUTPUT      output/cache/random_forest_inputs.RData   (inputs, clinical)
##              output/tables/SFig4_random_forest_input_summary.csv
##
##  METHODS     No statistics.  The matrices are the processed values exactly as
##              used for the submitted models (metagenomic profiles of all 42
##              stool samples; 31 fecal metabolomes).  They are stored "frozen"
##              because the metabolite matrix carries an assay-floor
##              subtraction (28 tiny negative residuals, minimum -5e-9) and 14
##              KO cells that differ from the current HUMAnN export - both are
##              retained so that the screening and models replay exactly.
##              To build the matrices from the raw tables instead, see the
##              comment block at the end of this script.
##
##  R PACKAGES  base R only
## =============================================================================

source("R/_common.R")
report_versions(character(0))


## ---- 1. Load the frozen matrices ----------------------------------------------

rf <- readRDS(input_file("random_forest", "rf_inputs_frozen.rds"))
inputs   <- rf$inputs
clinical <- rf$clinical
cat(rf$note, "\n")

stopifnot(identical(names(inputs), c("Strain", "Species", "Genus", "Family", "Phylum", "Pathway", "KO", "Metabolite")))
for (nm in names(inputs)) {
  x <- inputs[[nm]]
  stopifnot(
    is.numeric(x$values), all(is.finite(x$values)),
    !anyDuplicated(rownames(x$values)), !anyDuplicated(colnames(x$values)),
    identical(rownames(x$values), x$metadata$SampleID),
    all(x$metadata$Response %in% response_levels),
    all(x$metadata$Timepoint %in% c("Before", "Ongoing"))
  )
}
stopifnot(nrow(clinical) == 26, !anyDuplicated(clinical$SubjectID))


## ---- 2. Summary table -----------------------------------------------------------

summary_tbl <- data.frame(
  Input      = names(inputs),
  Samples    = vapply(inputs, function(x) nrow(x$values), integer(1)),
  Patients   = vapply(inputs, function(x) length(unique(x$metadata$SubjectID)), integer(1)),
  Baseline   = vapply(inputs, function(x) sum(x$metadata$Timepoint == "Before"), integer(1)),
  After_RT   = vapply(inputs, function(x) sum(x$metadata$Timepoint == "Ongoing"), integer(1)),
  Features   = vapply(inputs, function(x) ncol(x$values), integer(1)),
  row.names  = NULL
)
print(summary_tbl, row.names = FALSE)
save_table(summary_tbl, "SFig4_random_forest_input_summary")

save(inputs, clinical, file = cache_file("random_forest_inputs.RData"))
message("Done: 09a_SFig4_Random_forest_inputs.R")


## ---- How the matrices relate to the public tables -------------------------------
## Strain/Species/Genus/Family/Phylum  = input/microbiome/metaphlan_<rank>.csv
##      (relative abundance %, all 42 samples; features with all-zero columns
##      removed; SGB names as "Species|SGBxxxx").  These five match the CSVs exactly.
## Pathway   = input/microbiome/humann_pathway_abundance.tsv, unstratified rows,
##      UNMAPPED/UNINTEGRATED removed, no coverage mask.
## KO        = input/microbiome/humann_ko_cpm.tsv, unstratified rows (14 cells
##      differ from the frozen matrix; provenance not recoverable).
## Metabolite = input/metabolome/fecal_metabolites.csv minus one after-RT sample
##      (31 of 32 observations), assay floor subtracted per metabolite.
