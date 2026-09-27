## =============================================================================
##  18_Fig5_Four_omics_preparation.R
##  ---------------------------------------------------------------------------
##  PANELS      (no figure) - matched, preprocessed blocks for Figure 5A
##              (Procrustes / generalized Procrustes) and the MOFA scripts 20-23
##
##  PURPOSE     Bring the four data layers (species, KO, fecal metabolites, tumour
##              transcriptome) to a common sample key (SubjectID__Timepoint),
##              apply layer-specific filtering / transformation, and compute the
##              ordinations and pairwise sample overlaps.
##
##  INPUT       input/microbiome/metaphlan_species.csv, humann_ko_cpm.tsv
##              input/metabolome/fecal_metabolites.csv
##              input/host_rnaseq/tumor_rna_counts.csv  *** RESTRICTED ***
##              input/metadata/stool_samples.csv, rna_samples.csv, patients.csv
##
##  OUTPUT      output/cache/four_omics_inputs.RData  (list `coherence`)
##                $modalities[[layer]]  sample_meta, relative/table matrices, dist, pcoa_scores
##                $pairs[[pair]]        overlapping samples, pair-specific PCoA scores
##                $four_block           the 22 samples with all four layers
##              output/tables/Fig5A_layer_summary.csv, Fig5A_pair_summary.csv,
##              Fig5A_metabolite_QC.csv
##
##  METHODS
##   * Species and KO (42 metagenomes): features present in >= 20 % of samples;
##     re-normalised to relative abundance; Bray-Curtis dissimilarity (vegan);
##     the MOFA table is the Hellinger (square-root) transform, z-scored.
##     KOs annotated as ribosomal proteins and UNGROUPED/UNMAPPED are removed.
##   * Metabolites (32 samples): repeated-minimum values are treated as the
##     assay floor and replaced by floor/2; metabolites kept when >= 50 % of
##     samples are finite, >= 50 % lie above the repeated floor, >= 2 distinct
##     detected values and non-zero variance remain (25 of 40); z-scored;
##     Euclidean distance.
##   * Host (27 biopsies): genes with >= 10 counts in >= 20 % of biopsies,
##     blind DESeq2 VST; the highest-variance genes accounting for 50 % of the
##     summed gene-wise variance are kept (no further scaling); Euclidean distance.
##   * Ordination: principal coordinates (ape::pcoa), up to 5 axes with
##     positive eigenvalues, recomputed for each pairwise sample overlap.
##
##  R PACKAGES  vegan, ape, DESeq2, dplyr
## =============================================================================

source("R/_common.R")
suppressPackageStartupMessages({
  library(vegan)
  library(ape)
  library(DESeq2)
  library(dplyr)
})
report_versions(c("vegan", "ape", "DESeq2"))

count_file <- input_file("host_rnaseq", "tumor_rna_counts.csv")
if (!file.exists(count_file)) stop("input/host_rnaseq/tumor_rna_counts.csv (restricted) is required for Figure 5.", call. = FALSE)

## ---- 0. Settings ---------------------------------------------------------------------

prevalence_cut <- 0.20          # species and KO
met_min_finite <- 0.50; met_min_detected <- 0.50; met_floor_repeats <- 2
host_min_count <- 10; host_min_fraction <- 0.20; host_cum_variance <- 0.50
pcoa_axes <- 5
time_levels <- c("Before", "Ongoing")

sample_key <- function(subject, timepoint) paste(subject, timepoint, sep = "__")
pcoa_scores <- function(d, k = pcoa_axes) {
  p <- pcoa(d); pos <- which(p$values$Eigenvalues > 1e-8)
  s <- p$vectors[, pos[seq_len(min(k, length(pos)))], drop = FALSE]
  colnames(s) <- paste0("Axis", seq_len(ncol(s))); s
}


## ---- 1. Sample tables ----------------------------------------------------------------

patients <- read_patients()
stool <- read_stool_samples() %>% transmute(OmicSampleID = SampleID, SubjectID, Timepoint, Response)
metab_all <- read_metabolites()
met_meta <- metab_all %>% transmute(OmicSampleID = MetaboliteSampleID, SubjectID, Timepoint, Response)
rna <- read_csv_na(input_file("metadata", "rna_samples.csv")) %>% inner_join(patients, by = "SubjectID") %>%
  transmute(OmicSampleID = RNA_sample_id, SubjectID, Timepoint, Response)

order_meta <- function(m) {
  m <- m %>% arrange(SubjectID, match(Timepoint, time_levels)) %>% mutate(SampleID = sample_key(SubjectID, Timepoint))
  stopifnot(!anyDuplicated(m$SampleID)); m
}
stool <- order_meta(stool); met_meta <- order_meta(met_meta); rna <- order_meta(rna)


## ---- 2. Compositional blocks: species and KO -----------------------------------------

compositional_block <- function(raw, meta, name) {
  raw <- raw[meta$OmicSampleID, , drop = FALSE]; rownames(raw) <- meta$SampleID
  keep <- colMeans(raw > 0) >= prevalence_cut & colSums(raw) > 0
  filtered <- raw[, keep, drop = FALSE]
  relative <- filtered / rowSums(filtered)
  hellinger <- sqrt(relative)
  hellinger <- hellinger[, apply(hellinger, 2, var) > 0, drop = FALSE]
  d <- vegdist(relative[, colnames(hellinger)], method = "bray")
  cat(sprintf("%-10s %d samples x %d features (prevalence >= %.0f%%)\n", name, nrow(raw), ncol(hellinger), 100 * prevalence_cut))
  list(name = name, sample_meta = meta, raw = raw, relative = relative[, colnames(hellinger)], table = scale(hellinger),
       dist = d, pcoa_scores = pcoa_scores(d), n_features = ncol(hellinger))
}

species_raw <- read_metaphlan("species", samples = stool$OmicSampleID, drop_empty = TRUE)
ko_tab <- read.delim(input_file("microbiome", "humann_ko_cpm.tsv"), check.names = FALSE, stringsAsFactors = FALSE, quote = "")
ko_tab <- ko_tab[!grepl("^UNGROUPED|^NO_NAME|^UNMAPPED", ko_tab$KO) &
                 !grepl("large subunit ribosomal protein|small subunit ribosomal protein", ko_tab$KO_name, ignore.case = TRUE), ]
ko_raw <- t(as.matrix(ko_tab[, stool$OmicSampleID])); colnames(ko_raw) <- ko_tab$KO
ko_annotation <- data.frame(KO = ko_tab$KO, KO_name = sub(" \\[EC:.*\\]$", "", ko_tab$KO_name),
                            EC = sub(".*\\[EC:([^]]+)\\].*", "\\1", ko_tab$KO_name), stringsAsFactors = FALSE)
ko_annotation$EC[!grepl("\\[EC:", ko_tab$KO_name)] <- NA

species <- compositional_block(species_raw, stool, "species")
ko      <- compositional_block(ko_raw, stool, "ko")


## ---- 3. Metabolite block ---------------------------------------------------------------

met_raw <- as.matrix(metab_all[, metabolite_columns()]); rownames(met_raw) <- met_meta$SampleID[match(metab_all$MetaboliteSampleID, met_meta$OmicSampleID)]
met_raw <- met_raw[met_meta$SampleID, ]

## repeated-floor QC: the minimum value, when it occurs >= 2 times (or is 0), is
## the detection floor; floor values are replaced by floor/2 for the analysis
met_qc <- do.call(rbind, lapply(colnames(met_raw), function(m) {
  x <- met_raw[, m]; floor <- min(x); at_floor <- abs(x - floor) <= max(abs(floor) * 1e-8, 1e-12)
  repeated <- floor <= 1e-12 || sum(at_floor) >= met_floor_repeats
  detected <- if (repeated) !at_floor else rep(TRUE, length(x))
  replacement <- if (!repeated) NA else if (floor > 0) floor / 2 else if (any(x[detected] > 0)) min(x[detected & x > 0]) / 2 else 0
  x2 <- x; if (repeated) x2[at_floor] <- replacement
  data.frame(Metabolite = m, floor_value = floor, floor_repeated = repeated, floor_fraction = if (repeated) mean(at_floor) else 0,
             detected_fraction = mean(detected), n_unique_detected = length(unique(x[detected])), variance = var(x2),
             keep = mean(is.finite(x)) >= met_min_finite & mean(detected) >= met_min_detected & length(unique(x[detected])) >= 2 & var(x2) > 0)
}))
save_table(met_qc, "Fig5A_metabolite_QC")
met_matrix <- met_raw
for (m in colnames(met_raw)) if (met_qc$floor_repeated[met_qc$Metabolite == m]) {
  x <- met_raw[, m]; at_floor <- abs(x - min(x)) <= max(abs(min(x)) * 1e-8, 1e-12)
  met_matrix[at_floor, m] <- if (min(x) > 0) min(x) / 2 else if (any(x[!at_floor] > 0)) min(x[!at_floor & x > 0]) / 2 else 0
}
met_matrix <- met_matrix[, met_qc$Metabolite[met_qc$keep]]
cat(sprintf("metabolite %d samples x %d metabolites retained of %d\n", nrow(met_matrix), ncol(met_matrix), nrow(met_qc)))
met_table <- scale(met_matrix)
metabolite <- list(name = "metabolite", sample_meta = met_meta, raw = met_raw, table = met_table,
                   dist = dist(met_table), pcoa_scores = pcoa_scores(dist(met_table)), n_features = ncol(met_table))


## ---- 4. Host block ------------------------------------------------------------------------

counts_tab <- read.csv(count_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_name <- with(counts_tab, ifelse(!is.na(Gene_Symbol) & Gene_Symbol != "", Gene_Symbol, ifelse(!is.na(Gene_ID) & Gene_ID != "", Gene_ID, Transcript_ID)))
counts <- rowsum(as.matrix(counts_tab[, rna$OmicSampleID]), gene_name); counts <- round(counts); mode(counts) <- "integer"
counts <- counts[rowSums(counts >= host_min_count) >= ceiling(host_min_fraction * ncol(counts)), ]
col_data <- data.frame(intercept = rep(1L, ncol(counts)), row.names = colnames(counts))   # no design variable: blind VST
dds <- DESeqDataSetFromMatrix(counts, col_data, design = ~ 1)
vst_all <- t(assay(varianceStabilizingTransformation(dds, blind = TRUE)))        # biopsies x genes
rownames(vst_all) <- rna$SampleID

gene_var <- sort(apply(vst_all, 2, var), decreasing = TRUE); gene_var <- gene_var[gene_var > 0]
cum_frac <- cumsum(gene_var) / sum(gene_var)
host_top_n <- which(cum_frac >= host_cum_variance)[1]                            # <-- 50 % cumulative variance rule
host_table <- vst_all[, names(gene_var)[seq_len(host_top_n)]]
cat(sprintf("host       %d biopsies x %d genes (of %d after count QC; %.0f%% cumulative variance)\n",
            nrow(host_table), ncol(host_table), length(gene_var), 100 * host_cum_variance))
host <- list(name = "host", sample_meta = rna, table = host_table, dist = dist(host_table),
             pcoa_scores = pcoa_scores(dist(host_table)), n_features = ncol(host_table),
             vst_qc = vst_all[, names(gene_var)],                                          # all count-QC genes (used by MOFA, script 20)
             gene_variance = data.frame(Gene = names(gene_var), variance = gene_var, cumulative_fraction = cum_frac, row.names = NULL))

modalities <- list(species = species, ko = ko, metabolite = metabolite, host = host)


## ---- 5. Pairwise overlaps and four-block set -------------------------------------------------

pair_block <- function(x, y) {
  ids <- intersect(x$sample_meta$SampleID, y$sample_meta$SampleID)
  meta <- x$sample_meta %>% filter(SampleID %in% ids) %>% arrange(SubjectID, match(Timepoint, time_levels))
  ids <- meta$SampleID
  block_dist <- function(b) if (b$name %in% c("species", "ko")) vegdist(b$relative[ids, ], method = "bray") else dist(b$table[ids, ])
  dx <- block_dist(x); dy <- block_dist(y)
  out <- list(pair_name = paste(x$name, y$name, sep = "_"), block_names = c(x$name, y$name), sample_meta = meta,
              n = length(ids), n_baseline = sum(meta$Timepoint == "Before"), n_subjects = n_distinct(meta$SubjectID))
  out[[paste0(x$name, "_pcoa_scores")]] <- pcoa_scores(dx); out[[paste0(y$name, "_pcoa_scores")]] <- pcoa_scores(dy)
  out[[paste0(x$name, "_table")]] <- x$table[ids, ]; out[[paste0(y$name, "_table")]] <- y$table[ids, ]
  out
}
pair_names <- combn(names(modalities), 2, paste, collapse = "_")
pairs <- lapply(combn(names(modalities), 2, simplify = FALSE), function(nm) pair_block(modalities[[nm[1]]], modalities[[nm[2]]]))
names(pairs) <- pair_names
pair_summary <- bind_rows(lapply(pairs, function(p) data.frame(pair = p$pair_name, n_samples = p$n, n_baseline = p$n_baseline, n_subjects = p$n_subjects)))
print(pair_summary, row.names = FALSE)
save_table(pair_summary, "Fig5A_pair_summary")

common <- Reduce(intersect, lapply(modalities, function(b) b$sample_meta$SampleID))
four_meta <- species$sample_meta %>% filter(SampleID %in% common)
four_block <- list(block_names = names(modalities), sample_meta = four_meta,
                   tables = lapply(modalities, function(b) b$table[four_meta$SampleID, ]))
cat("Samples with all four layers:", nrow(four_meta), "(", sum(four_meta$Timepoint == "Before"), "baseline )\n")

layer_summary <- bind_rows(lapply(modalities, function(b) data.frame(layer = b$name, n_samples = nrow(b$sample_meta), n_features = b$n_features)))
save_table(layer_summary, "Fig5A_layer_summary")

coherence <- list(modalities = modalities, pairs = pairs, four_block = four_block, ko_annotation = ko_annotation,
                  settings = list(prevalence_cut = prevalence_cut, met_min_finite = met_min_finite, met_min_detected = met_min_detected,
                                  host_min_count = host_min_count, host_min_fraction = host_min_fraction, host_cum_variance = host_cum_variance,
                                  host_top_n = host_top_n, pcoa_axes = pcoa_axes))
save(coherence, file = cache_file("four_omics_inputs.RData"))
message("Done: 18_Fig5_Four_omics_preparation.R")
