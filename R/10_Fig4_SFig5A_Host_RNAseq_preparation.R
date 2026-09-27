## =============================================================================
##  10_Fig4_SFig5A_Host_RNAseq_preparation.R
##  ---------------------------------------------------------------------------
##  PANELS      Supplementary Figure 5A  volcano plot of tumour gene expression,
##                                       pCR versus non-pCR (baseline biopsies)
##              (+ the expression objects used by scripts 11-17 and 18-23)
##
##  PURPOSE     Differential expression of baseline tumour biopsies between
##              response groups, and a variance-stabilised expression matrix
##              for the correlation / integration analyses.
##
##  INPUT       input/host_rnaseq/tumor_rna_counts.csv   *** RESTRICTED ***
##                  gene-level read counts (Gene_ID, Gene_Symbol, Transcript_ID,
##                  gene_biotype, one column per RNA sample); not part of the
##                  public deposit - scripts 10 and 13-17 need it, scripts 11/12
##                  can run from the frozen GSEA tables instead.
##              input/metadata/rna_samples.csv        27 biopsies -> patient / visit
##              input/metadata/patients.csv
##              input/metabolome/fecal_metabolites.csv (to define the RNA-metabolome
##                  matched subset used by scripts 14-17)
##
##  OUTPUT      output/figures/SFig5A_host_volcano.{svg,png}
##              output/tables/SFig5A_host_DESeq2_results.csv
##              output/cache/host_deseq2.RData
##                  dds        DESeqDataSet (20 baseline biopsies)
##                  deg        results table (log2FC pCR / non-pCR, stat, P, FDR, Wilcoxon P)
##                  vst        VST expression matrix (genes x 20 biopsies)
##                  rna_meta   biopsy table in matrix column order
##                  matched    the 17 biopsies with a baseline fecal metabolome
##                  met_matrix metabolite concentrations of those 17 patients
##
##  METHODS     * Genes: rows collapsed to one per gene symbol (Gene_Symbol, else
##                Gene_ID, else Transcript_ID), counts summed; genes with >= 10
##                reads in >= 3 biopsies are kept (22,077).
##              * DESeq2 1.42.1: negative-binomial GLM, design ~ Response
##                (non-pCR reference), Wald test; log2FC > 0 = higher in pCR;
##                FDR = Benjamini-Hochberg (alpha = 0.1 for independent filtering).
##              * VST = DESeq2 variance-stabilising transformation (blind = FALSE)
##                for correlations; a per-gene Wilcoxon rank-sum test on VST
##                values is added as a rank-based check.
##              * Volcano (S5A): colour = direction for nominal P < 0.05; the 20
##                genes with the lowest FDR are labelled.
##
##  R PACKAGES  DESeq2 1.42.1, SummarizedExperiment, dplyr, tidyr, ggplot2, ggrepel
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
})
report_versions(c("DESeq2", "ggplot2", "ggrepel"))

count_file <- input_file("host_rnaseq", "tumor_rna_counts.csv")
if (!file.exists(count_file)) {
  stop("input/host_rnaseq/tumor_rna_counts.csv (restricted host RNA-seq counts) is not available.\n",
       "Scripts 11 and 12 can still be run from the frozen GSEA tables; scripts 13-17 need this file.", call. = FALSE)
}


## ---- 1. Baseline tumour biopsies ------------------------------------------------

rna_samples <- read_csv_na(input_file("metadata", "rna_samples.csv"))
rna_meta <- rna_samples %>%
  dplyr::filter(Timepoint == "Before") %>%
  inner_join(read_patients(), by = "SubjectID") %>%
  arrange(EnrolmentOrder)
rna_meta$Response <- factor(rna_meta$Response, levels = response_levels)
rna_meta$Group <- factor(ifelse(rna_meta$Response == "pCR", "pCR", "non_pCR"), levels = c("non_pCR", "pCR"))  # DESeq2 reference = non_pCR
rownames(rna_meta) <- rna_meta$RNA_sample_id
cat("Baseline tumour biopsies:", nrow(rna_meta), "(pCR", sum(rna_meta$Response == "pCR"), "/ non-pCR",
    sum(rna_meta$Response == "non-pCR"), ")\n")

## The subset with a baseline fecal metabolome (used by scripts 14-17)
metab <- read_metabolites(timepoint = "Before")
matched <- rna_meta %>% dplyr::filter(SubjectID %in% metab$SubjectID)
met_matrix <- as.matrix(metab[match(matched$SubjectID, metab$SubjectID), metabolite_columns()])
rownames(met_matrix) <- matched$RNA_sample_id
cat("Biopsies with a matched baseline metabolome:", nrow(matched), "\n")


## ---- 2. Count matrix (one row per gene) -----------------------------------------

rna <- read.csv(count_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(rna_meta$RNA_sample_id %in% names(rna)))

gene_name <- with(rna, ifelse(!is.na(Gene_Symbol) & Gene_Symbol != "", Gene_Symbol,
                              ifelse(!is.na(Gene_ID) & Gene_ID != "", Gene_ID, Transcript_ID)))
counts <- rowsum(as.matrix(rna[, rna_meta$RNA_sample_id]), gene_name)        # sum rows sharing a symbol
counts <- round(counts); mode(counts) <- "integer"
cat("Genes in the count table:", nrow(counts), "\n")


## ---- 3. DESeq2 -----------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = counts, colData = rna_meta, design = ~ Group)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]                          # <-- low-count gene filter
cat("Genes retained after filtering:", nrow(dds), "\n")

dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("Group", "pCR", "non_pCR"), alpha = 0.10)

vst <- assay(vst(dds, blind = FALSE))                                   # genes x biopsies
stopifnot(identical(colnames(vst), rna_meta$RNA_sample_id))

## rank-based check on the VST values
wilcox_vst <- apply(vst, 1, function(x) wilcox_p(x, rna_meta$Response, exact = FALSE))

deg <- data.frame(Gene = rownames(res), as.data.frame(res), row.names = NULL) %>%
  dplyr::rename(P = pvalue, FDR = padj) %>%                              # dplyr:: avoids Bioconductor name clashes
  mutate(Wilcoxon_P = wilcox_vst[Gene], Wilcoxon_FDR = p.adjust(Wilcoxon_P, "BH"),
         Direction = ifelse(log2FoldChange > 0, "pCR high", "non-pCR high")) %>%
  arrange(P)
cat(sprintf("Nominal P < 0.05: %d pCR-high / %d non-pCR-high;  FDR < 0.05: %d / %d\n",
            sum(deg$P < 0.05 & deg$log2FoldChange > 0, na.rm = TRUE), sum(deg$P < 0.05 & deg$log2FoldChange < 0, na.rm = TRUE),
            sum(deg$FDR < 0.05 & deg$log2FoldChange > 0, na.rm = TRUE), sum(deg$FDR < 0.05 & deg$log2FoldChange < 0, na.rm = TRUE)))
save_table(deg, "SFig5A_host_DESeq2_results")

save(dds, deg, vst, rna_meta, matched, met_matrix, file = cache_file("host_deseq2.RData"))


## ---- 4. Supplementary Figure 5A: volcano plot -----------------------------------

deg$Class <- factor(case_when(is.na(deg$P) | deg$P >= 0.05 ~ "Not significant",
                              deg$log2FoldChange > 0 ~ "pCR high", TRUE ~ "non-pCR high"),
                    levels = c("Not significant", "pCR high", "non-pCR high"))
label_genes <- deg %>% dplyr::filter(!is.na(FDR)) %>% arrange(FDR, desc(abs(log2FoldChange))) %>% slice_head(n = 20)   # <-- labelled genes

p_volcano <- ggplot(dplyr::filter(deg, !is.na(P)), aes(log2FoldChange, -log10(P), colour = Class)) +
  geom_point(size = 0.9, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_text_repel(data = label_genes, aes(label = Gene), size = 2.4, max.overlaps = 100, seed = 1,
                  box.padding = 0.3, segment.colour = "grey60", show.legend = FALSE) +
  scale_colour_manual(values = c("Not significant" = "#C8C8C8", "pCR high" = response_colors[["pCR"]],
                                 "non-pCR high" = response_colors[["non-pCR"]]), name = NULL) +
  labs(x = expression(Log[2] * " fold change (pCR / non-pCR)"), y = expression(-Log[10] * "(" * italic(p) * "-value)")) +
  theme_tnt(base_size = 10, legend = "bottom") +
  theme(axis.line = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6))
save_panel(p_volcano, "SFig5A_host_volcano", width = 4.6, height = 5)

message("Done: 10_Fig4_SFig5A_Host_RNAseq_preparation.R")
