## =============================================================================
##  14_Fig4DE_SFig5B_Cross_omics_preparation.R
##  ---------------------------------------------------------------------------
##  PANELS      (no figure) - quality-controlled metabolite and host-gene tables
##              for Figure 4D, Figure 4E and Supplementary Figure 5B
##
##  PURPOSE     Define the 17 patients with both a baseline fecal metabolome and
##              a baseline tumour biopsy, remove metabolites that are essentially
##              constant / at the assay floor in this subset, and pre-select the
##              host genes that enter the metabolite-host correlation analyses.
##
##  INPUT       output/cache/host_deseq2.RData   vst, deg, matched, met_matrix (script 10)
##
##  OUTPUT      output/cache/host_metabolite_inputs.RData
##                  meta          17 patients (RNA_sample_id, SubjectID, Response ...)
##                  met_raw       metabolite concentrations (umol/g) after QC
##                  met_log       log10(concentration + 1e-6)
##                  host_expr     VST expression of the selected genes (patients x genes)
##                  host_deg      DESeq2 / Wilcoxon statistics of all genes
##                  immune_genes  curated immune gene list
##              output/tables/Fig4D_metabolite_QC.csv
##
##  METHODS
##   * Metabolite QC (17 patients): >= 8 finite values, >= 4 distinct values,
##     SD and IQR > 0, < 65 % of samples at the exact minimum (assay floor) and
##     < 50 % within a "basal window" of 5 % of the IQR above the minimum.
##   * Host genes: (i) genes with DESeq2 P < 0.1 or VST Wilcoxon P < 0.1 in the
##     17 patients that pass a general QC (>= 4 distinct values, SD and IQR >
##     0, < 80 % at the floor), plus (ii) a curated immune gene list (T cell,
##     NK, exhaustion, MDSC, neutrophil, antigen presentation, TLS/B cell ...)
##     under a relaxed QC, regardless of differential expression.
##
##  R PACKAGES  base R
## =============================================================================

source("R/_common.R")
report_versions(character(0))

load(cache_file("host_deseq2.RData"))                      # vst, deg, matched, met_matrix
meta <- matched
stopifnot(identical(rownames(met_matrix), meta$RNA_sample_id), all(meta$RNA_sample_id %in% colnames(vst)))
cat("Patients with baseline biopsy + fecal metabolome:", nrow(meta), "\n")

host_expr_all <- t(vst[, meta$RNA_sample_id])              # patients x genes


## ---- 1. Metabolite quality control ----------------------------------------------

met_qc <- do.call(rbind, lapply(colnames(met_matrix), function(m) {
  x <- met_matrix[, m]; x <- x[is.finite(x)]
  data.frame(Metabolite = m, n_finite = length(x), n_distinct = length(unique(round(x, 8))),
             SD = sd(x), IQR = IQR(x),
             floor_fraction = mean(round(x, 8) == round(min(x), 8)),                 # share at the exact minimum
             basal_fraction = mean(x <= min(x) + max(1e-8, 0.05 * IQR(x))))          # share within the basal window
}))
met_qc$keep <- with(met_qc, n_finite >= 8 & n_distinct >= 4 & SD > 0 & IQR > 0 & floor_fraction < 0.65 & basal_fraction < 0.5)   # <-- QC rule
save_table(met_qc, "Fig4D_metabolite_QC")
met_raw <- met_matrix[, met_qc$Metabolite[met_qc$keep], drop = FALSE]
met_log <- log10(met_raw + 1e-6)
cat("Metabolites passing QC:", ncol(met_raw), "of", ncol(met_matrix), "\n")


## ---- 2. Host gene selection ------------------------------------------------------

## Curated immune / tumour-microenvironment genes (kept regardless of DE evidence)
immune_genes <- unique(c(
  "CD3D", "CD3E", "TRAC", "CD8A", "CD8B", "PRF1", "GZMA", "GZMB", "GZMH", "NKG7", "CTSW", "CCL5", "IFNG", "TNF", "IL2", "LTA",
  "CCL3", "CCL4", "XCL1", "XCL2", "CD69", "CD38", "TNFRSF9", "IL2RA", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "ENTPD1", "TOX", "CXCL13", "LAYN",
  "TCF7", "SLAMF6", "CXCR5", "IL7R", "CCR7", "LEF1", "ITGAE", "CXCR6", "ITGA1", "ZNF683", "RUNX3", "MKI67", "TOP2A", "STMN1", "TYMS", "PCNA",
  "CDK1", "CCNB1", "CCNB2", "FASLG", "TNFSF10", "FAS", "TNFRSF10A", "TNFRSF10B", "TNFRSF1A", "CD274", "PDCD1LG2", "CD80", "CD86", "PVR",
  "NECTIN2", "LGALS9", "FGL1", "STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "CXCR3", "FOXP3", "CTLA4", "IKZF2", "TNFRSF18", "CCR8", "ICOS",
  "IL10", "TGFB1", "EBI3", "IL12A", "NT5E", "LRRC32", "FGL2", "OLR1", "S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "ARG1", "CEACAM8",
  "CD14", "CCR2", "VCAN", "FCN1", "LILRB1", "IL1B", "IDO1", "PTGS2", "CYBB", "NCF1", "NCF2", "NCF4", "RORC", "IL17A", "IL17F", "CCR6",
  "KLRB1", "IL23R", "RORA", "CCL20", "IL1R1", "CSF2", "TBX21", "BHLHE40", "IL17RA", "IL17RC", "TRAF3IP2", "NFKBIZ", "CXCL1", "CXCL2",
  "CXCL5", "CXCL6", "CXCL8", "CSF3", "KLRD1", "NCR1", "KLRK1", "FCGR3A", "TYROBP", "FCER1G", "GNLY", "CLEC9A", "XCR1", "BATF3", "IRF8",
  "WDFY4", "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5", "PSMB8", "PSMB9", "MS4A1", "CD79A", "CD79B", "CD74",
  "CCL19", "CCL21", "MZB1", "JCHAIN", "IGKC", "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "SPP1", "CD163", "MRC1", "MARCO", "CSF1R",
  "GATA3", "IL4", "IL5", "IL13", "PTGDR2", "CEBPB", "STAT3", "KLRC1", "KLRC2", "KLRF1", "FGFBP2", "CD244", "CD160", "TOX2", "NR4A1",
  "NR4A2", "NR4A3", "EOMES", "PRDM1", "IL12RB2", "IL18R1", "IL4R", "STAT6", "IL21", "IL22", "IL26", "AHR", "FPR1", "FPR2", "CXCR1",
  "S100A12", "MMP8", "MMP9", "LTF", "CAMP", "LCN2", "OLFM4", "MNDA", "SELL", "SLC7A2", "CADM1", "THBD", "DNASE1L3", "SNX22", "CPVL",
  "CLNK", "BCL6", "IL21R", "CD37", "BANK1"
))

host_qc <- data.frame(
  Gene = colnames(host_expr_all),
  N = colSums(is.finite(host_expr_all)),
  Distinct = apply(host_expr_all, 2, function(x) length(unique(round(x, 6)))),
  SD = apply(host_expr_all, 2, sd), IQR = apply(host_expr_all, 2, IQR),
  Floor_fraction = apply(host_expr_all, 2, function(x) mean(round(x, 6) == round(min(x), 6))),
  row.names = NULL
)
host_qc$General <- with(host_qc, N >= 8 & Distinct >= 4 & SD > 0 & IQR > 0 & Floor_fraction < 0.8)    # <-- QC for DE-eligible genes
host_qc$Curated <- with(host_qc, N >= 8 & Distinct >= 3 & SD > 0 & Floor_fraction < 0.95)             # <-- relaxed QC for immune genes

host_deg <- data.frame(Gene = deg$Gene, LFC = deg$log2FoldChange, DESeq2_P = deg$P, DESeq2_FDR = deg$FDR)
host_deg$Wilcoxon_P <- apply(host_expr_all[, host_deg$Gene], 2, function(x) wilcox_p(x, meta$Response, exact = FALSE))   # 17 patients
host_deg$Eligible <- with(host_deg, (!is.na(DESeq2_P) & DESeq2_P < 0.1) | (!is.na(Wilcoxon_P) & Wilcoxon_P < 0.1))    # <-- DE evidence

qc <- host_qc[match(host_deg$Gene, host_qc$Gene), ]
selected <- (host_deg$Eligible & qc$General) | (host_deg$Gene %in% immune_genes & qc$Curated)
host_expr <- host_expr_all[, host_deg$Gene[selected], drop = FALSE]
cat("Host genes selected:", ncol(host_expr), "(", sum(host_deg$Eligible & qc$General), "DE-eligible +",
    sum(host_deg$Gene %in% immune_genes & qc$Curated & !(host_deg$Eligible & qc$General)), "curated immune )\n")

save(meta, met_raw, met_log, host_expr, host_deg, host_qc, met_qc, immune_genes,
     file = cache_file("host_metabolite_inputs.RData"))
message("Done: 14_Fig4DE_SFig5B_Cross_omics_preparation.R")
