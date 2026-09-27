## =============================================================================
##  R/_common.R  --  shared settings for every analysis script
##  ---------------------------------------------------------------------------
##  This file is sourced at the top of each numbered script:
##
##      source("R/_common.R")
##
##  It (1) checks that R is running from the package root, (2) defines the
##  folder layout, (3) defines the colour palettes and ggplot2 theme used in
##  every figure, and (4) provides a handful of small helper functions
##  (metadata loading, Wilcoxon tests, P-value formatting, figure export).
##
##  Nothing here performs an analysis.  Change colours / theme / output folder
##  here and every panel follows.
## =============================================================================

## ---- 1. Package root and folder layout ---------------------------------------
## The scripts are written to be run with the package folder as the working
## directory (open LARC_TNT_response.Rproj in RStudio, or setwd() to it).

if (!file.exists("LARC_TNT_response.Rproj") || !dir.exists("input")) {
  stop(
    "Please run the scripts from the package root (the folder that contains\n",
    "LARC_TNT_response.Rproj and input/).  In RStudio: File > Open Project...\n",
    "or setwd('<path>/LARC_TNT_response') first.", call. = FALSE
  )
}

paths <- list(
  root    = normalizePath(getwd(), winslash = "/"),
  input   = "input",                 # curated, ready-to-use input files (see input/README.md)
  figures = "output/figures",        # one SVG + one PNG per figure panel
  tables  = "output/tables",         # numeric results behind each panel (CSV)
  cache   = "output/cache"           # intermediate R objects passed between scripts
)
for (p in paths[c("figures", "tables", "cache")]) dir.create(p, recursive = TRUE, showWarnings = FALSE)

input_file <- function(...) file.path(paths$input, ...)     # e.g. input_file("metadata", "patients.csv")
cache_file <- function(...) file.path(paths$cache, ...)     # e.g. cache_file("host_deseq2.RData")

## ---- 2. Colours ---------------------------------------------------------------
## Response-group colours are taken from the submitted vector figures
## (pCR = teal #4FAE9A, non-pCR = salmon #DE7872).

response_colors <- c("pCR" = "#4FAE9A", "non-pCR" = "#DE7872")
response_levels <- c("pCR", "non-pCR")

## TRG score tiles (0 = white ... 3 = dark brown), used in heatmap annotations
trg_score_colors <- c("0" = "#FFFFFF", "1" = "#D9BF8C", "2" = "#B98341", "3" = "#7A5C3A")

## Family colours (Figure 1D, Supplementary Figure 1C)
family_colors <- c(
  "Others"                     = "#BEBEBE",
  "Lachnospiraceae"            = "#FF6B6B",
  "Oscillospiraceae"           = "#F7D9BC",
  "Bifidobacteriaceae"         = "#E0B0FF",
  "Bacteroidaceae"             = "#F88379",
  "Erysipelotrichaceae"        = "#FBEC5D",
  "Clostridiaceae"             = "#B5C7EB",
  "Prevotellaceae"             = "#00BB77",
  "Streptococcaceae"           = "#E89EB8",
  "Eubacteriales_unclassified" = "#6395EE",
  "Coriobacteriaceae"          = "#FFFF99",
  "Enterobacteriaceae"         = "#B5E48C",   # appear only in the per-group top-10 union (S1C)
  "Lactobacillaceae"           = "#F4E285",
  "Selenomonadaceae"           = "#C4A3FF"
)

## Curated microbial-group tile colours (Figure 1F/1G)
group_tile_colors <- c(
  "Direct butyrate producer" = "#8BCF8A",
  "Butyrate cross-feeder"    = "#C5E58A",
  "ICB-favorable bacteria"   = "#8FADE8",
  "CRC signatures"           = "#F28C8C",
  "Unassigned"               = "#E5E5E5"
)

## Diverging palettes used across heatmaps and correlation plots
zscore_colors      <- c("#4575B4", "#F7F7F7", "#D73027")   # blue - white - red (row z-scores)
correlation_colors <- c("#2166AC", "#FFFFFF", "#B2182B")   # blue - white - red (Spearman rho)
log2fc_colors      <- c("#DE7872", "#FFFFFF", "#4FAE9A")   # non-pCR - white - pCR
p_value_colors     <- c("#FFFFFF", "#2166AC")              # white - blue (-log10 P)

## ---- 3. ggplot2 theme ---------------------------------------------------------

suppressPackageStartupMessages(library(ggplot2))

## When several scripts are sourced in ONE R session, packages attached by an
## earlier script can mask the dplyr verbs used later (mixOmics attaches MASS,
## whose select() hides dplyr::select()).  Detaching MASS and re-attaching
## dplyr / tidyr on top of the search path (each script calls library(dplyr)
## after this file) keeps every script self-consistent.  Running the scripts in
## fresh sessions (run_all.R does) makes this a no-op.
for (p in c("MASS", "tidyr", "dplyr")) {
  if (paste0("package:", p) %in% search()) suppressWarnings(detach(paste0("package:", p), character.only = TRUE, force = TRUE))
}

## A publication theme close to ggpubr::theme_pubr(): black axis lines, no grid,
## legend on top.  `base_size` controls every font size of a panel at once.
theme_tnt <- function(base_size = 11, legend = "top") {
  theme_classic(base_size = base_size) +
    theme(
      axis.line        = element_line(colour = "black", linewidth = 0.5),
      axis.ticks       = element_line(colour = "black", linewidth = 0.5),
      axis.text        = element_text(colour = "black"),
      legend.position  = legend,
      legend.key       = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(hjust = 0.5)
    )
}

theme_set(theme_tnt())

## ---- 4. Metadata helpers ------------------------------------------------------

## read.csv with empty cells treated as missing (also for character columns)
read_csv_na <- function(file, ...) {
  read.csv(file, stringsAsFactors = FALSE, na.strings = c("", "NA"), ...)
}

## Patients (26) with the response label used everywhere: pCR / non-pCR.
## Rows are returned in EnrolmentOrder (an anonymous 1..26 index that
## reproduces the patient order of the original analyses).
##
## clinical = TRUE adds patient covariates from whichever file is present:
##   input/metadata/patients_clinical.csv    full clinical table (RESTRICTED,
##                                           available on request; 26 variables)
##   input/metadata/patients_covariates.csv  public minimum: Sex, Age, BMI - the
##                                           only covariates the adjusted
##                                           correlations (S5B, S6, 5F) use
## has_covariates() tells a script whether the columns it needs are available.
read_patients <- function(clinical = FALSE) {
  x <- read_csv_na(input_file("metadata", "patients.csv"))
  x$Response <- factor(x$Response, levels = response_levels)
  if (clinical) {
    f <- input_file("metadata", c("patients_clinical.csv", "patients_covariates.csv"))
    f <- f[file.exists(f)][1]
    if (!is.na(f)) x <- merge(x, read_csv_na(f), by = "SubjectID", all.x = TRUE, sort = FALSE)   # empty cells -> NA
  }
  x <- x[order(x$EnrolmentOrder), ]
  rownames(x) <- NULL
  x
}
has_covariates <- function(cols) all(cols %in% names(read_patients(clinical = TRUE)))

## Stool samples (42) joined to the patient table.
## Use `timepoint = "Before"` for the 26 baseline samples, "Ongoing" for the
## 16 after-RT samples, or NULL for all 42.
## ORDER MATTERS: rows come back in AnalysisOrder, the sample order of the
## original analysis.  Permutation tests (PERMANOVA, envfit) and bootstraps
## draw random numbers per row, so a different row order gives slightly
## different P values / intervals even with the same seed.
read_stool_samples <- function(timepoint = NULL, clinical = FALSE) {
  s <- read_csv_na(input_file("metadata", "stool_samples.csv"))
  s <- merge(s, read_patients(clinical = clinical), by = "SubjectID", sort = FALSE)
  if (!is.null(timepoint)) s <- s[s$Timepoint %in% timepoint, ]
  s <- s[order(s$AnalysisOrder), ]
  rownames(s) <- NULL
  s
}

## MetaPhlAn table (input/microbiome/metaphlan_<rank>.csv) -> samples x features
## matrix of relative abundance (%), feature names cleaned to the last rank.
read_metaphlan <- function(rank = c("species", "sgb", "genus", "family", "phylum"),
                           samples = NULL, drop_empty = TRUE) {
  rank <- match.arg(rank)
  x <- read.csv(input_file("microbiome", paste0("metaphlan_", rank, ".csv")),
                check.names = FALSE, stringsAsFactors = FALSE)
  feature <- switch(rank,
    sgb     = gsub("t__", "", sub(".*s__", "", x$clade_name)),   # "Species|SGBxxxx"
    species = sub(".*s__", "", x$clade_name),
    genus   = sub(".*g__", "", x$clade_name),
    family  = sub(".*f__", "", x$clade_name),
    phylum  = sub(".*p__", "", x$clade_name)
  )
  stopifnot(!anyDuplicated(feature))
  if (is.null(samples)) samples <- setdiff(names(x), "clade_name")
  m <- t(as.matrix(x[, samples, drop = FALSE]))
  colnames(m) <- feature
  if (drop_empty) m <- m[, colSums(m) > 0, drop = FALSE]
  m
}

## Generic reader for the HUMAnN / CAZy tables in input/microbiome/ (features in
## rows, samples in columns).  Returns a features x samples numeric matrix with
## the identifier column as row names; extra annotation columns are dropped.
##   read_feature_table("humann_ko_cpm.tsv", "KO", samples)
read_feature_table <- function(file, id_col, samples = NULL) {
  x <- read.delim(input_file("microbiome", file), check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  if (is.null(samples)) samples <- grep("^Sample_", names(x), value = TRUE)
  stopifnot(all(samples %in% names(x)), !anyDuplicated(x[[id_col]]))
  m <- as.matrix(x[, samples, drop = FALSE])
  rownames(m) <- x[[id_col]]
  m
}

## Fecal metabolite concentrations (umol/g), one row per metabolome sample.
## Rows are returned in AnalysisOrder (= the laboratory table order used by the
## original ordination code); see the note above read_stool_samples().
read_metabolites <- function(timepoint = NULL) {
  x <- read.csv(input_file("metabolome", "fecal_metabolites.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  x <- merge(x, read_patients(), by = "SubjectID", sort = FALSE)
  if (!is.null(timepoint)) x <- x[x$Timepoint %in% timepoint, ]
  x <- x[order(x$AnalysisOrder), ]
  rownames(x) <- NULL
  x
}

## The 40 metabolite column names, in assay order
metabolite_columns <- function() {
  x <- read.csv(input_file("metabolome", "fecal_metabolites.csv"), check.names = FALSE, nrows = 1)
  setdiff(names(x), c("AnalysisOrder", "MetaboliteSampleID", "SubjectID", "Timepoint", "SampleID"))
}

## ---- 4b. Taxon-name formatting ------------------------------------------------

## "Ruminococcus_sp_AF13_28|SGB4834" -> "Ruminococcus sp. AF13-28|SGB4834"
## (MetaPhlAn clade names use "_" both as a word separator and inside strain codes)
pretty_species <- function(x) {
  x <- gsub("_sp_", " sp. ", x)
  x <- gsub("_subsp_", " subsp. ", x)
  ## strain-code parts (upper-case letters/digits) keep a hyphen: AF13_28 -> AF13-28
  x <- gsub("(?<=[A-Z0-9])_(?![SG]GB)(?=[A-Z0-9])", "-", x, perl = TRUE)   # but "GGB9597_SGB15022" keeps a space
  gsub("_", " ", x)
}

## Which words of a species name are written upright (not italic)?
.upright_word <- function(words) {
  grepl("^(sp\\.|subsp\\.|bacterium|GGB|SGB|oral|taxon|group)", words) |
    grepl("[0-9]", words) | grepl("(aceae|ales|phylum)$", words) | grepl("^[A-Z]{2,}", words)
}

## plotmath label: Latin binomial in italics, everything else (sp., subsp.,
## strain codes, SGB/GGB identifiers, "bacterium") upright.  Returns a
## character vector of plotmath strings; use  labels = function(b) parse(text = lab[b])
## in a discrete ggplot2 scale, or parse(text = taxon_plotmath(x)) directly.
taxon_plotmath <- function(x) {
  vapply(pretty_species(x), function(s) {
    parts <- strsplit(s, "|", fixed = TRUE)[[1]]
    words <- strsplit(parts[1], " ", fixed = TRUE)[[1]]
    up <- .upright_word(words)
    run <- cumsum(c(TRUE, up[-1] != up[-length(up)]))         # runs of italic / upright words
    pieces <- vapply(split(seq_along(words), run), function(i) {
      txt <- paste(words[i], collapse = " ")
      if (up[i[1]]) sprintf('"%s"', txt) else sprintf('italic("%s")', txt)
    }, character(1))
    out <- paste(pieces, collapse = "~")
    if (length(parts) > 1) out <- paste0(out, '*"|', parts[2], '"')
    out
  }, character(1), USE.NAMES = FALSE)
}

## The same, as Markdown (for ggtext::element_markdown(), if that package is available)
taxon_markdown <- function(x) {
  vapply(pretty_species(x), function(s) {
    parts <- strsplit(s, "|", fixed = TRUE)[[1]]
    words <- strsplit(parts[1], " ", fixed = TRUE)[[1]]
    up <- .upright_word(words)
    words[!up] <- paste0("*", words[!up], "*")
    out <- gsub("\\* \\*", " ", paste(words, collapse = " "))
    if (length(parts) > 1) out <- paste0(out, "|", parts[2])
    out
  }, character(1), USE.NAMES = FALSE)
}

## ---- 5. Statistics helpers ----------------------------------------------------

## Two-sided Wilcoxon rank-sum P value for a numeric vector split by group.
## `exact = NULL` (the default) lets R decide, exactly as the original analyses did;
## pass exact = FALSE to force the normal approximation.
wilcox_p <- function(x, group, exact = NULL) {
  g <- as.character(group)
  suppressWarnings(wilcox.test(x[g == response_levels[1]], x[g == response_levels[2]], exact = exact)$p.value)
}

## Significance stars
p_stars <- function(p, cut = c(0.001, 0.01, 0.05), symbols = c("***", "**", "*")) {
  out <- rep("", length(p))
  for (i in seq_along(cut)) out[!is.na(p) & p < cut[i] & out == ""] <- symbols[i]
  out
}

## "P = 0.047" / "P < 0.001" style labels
p_label <- function(p, digits = 3, prefix = "P") {
  ifelse(p < 0.001, paste0(prefix, " < 0.001"), paste0(prefix, " = ", formatC(p, format = "f", digits = digits)))
}

## Hedges' g (bias-corrected standardised mean difference), positive = higher in x
hedges_g <- function(x, y) {
  n1 <- length(x); n2 <- length(y); df <- n1 + n2 - 2
  s <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / df)
  if (!is.finite(s) || s <= 0) return(NA_real_)
  ((mean(x) - mean(y)) / s) * (1 - 3 / (4 * df - 1))
}

## ---- 5b. Shared ComplexHeatmap pieces -----------------------------------------

## Legend font sizes used by every ComplexHeatmap panel (call after library(ComplexHeatmap))
set_heatmap_options <- function() {
  ComplexHeatmap::ht_opt(
    legend_title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    legend_labels_gp = grid::gpar(fontsize = 8),
    legend_grid_height = grid::unit(3.5, "mm"), legend_grid_width = grid::unit(3.5, "mm"),
    heatmap_row_names_gp = grid::gpar(fontsize = 8), message = FALSE
  )
}

## Top annotation used by every sample heatmap: TRG score tiles (0-3) and the
## response group.  `meta` must have columns TRG_score and Response, in the
## column order of the heatmap matrix.
response_top_annotation <- function(meta) {
  ComplexHeatmap::HeatmapAnnotation(
    `TRG score` = factor(meta$TRG_score, levels = 0:3), TRG = meta$Response,
    col = list(`TRG score` = trg_score_colors, TRG = response_colors),
    gp = grid::gpar(col = "grey40", lwd = 0.5), simple_anno_size = grid::unit(3.5, "mm"),
    annotation_name_gp = grid::gpar(fontsize = 9), annotation_name_side = "right",
    annotation_legend_param = list(`TRG score` = list(nrow = 2), TRG = list(nrow = 2))
  )
}

## Sample order for heatmaps: pCR block first, then non-pCR; within each block
## the samples are ordered by hierarchical clustering (Euclidean / complete) of
## the (z-scored) matrix `z` (features x samples).
order_samples_within_response <- function(z, meta, method = "complete") {
  unlist(lapply(response_levels, function(r) {
    ids <- meta$SampleID[meta$Response == r]
    if (length(ids) < 3) return(ids)
    ids[hclust(dist(t(z[, ids, drop = FALSE])), method = method)$order]
  }))
}

## ---- 6. Figure export ---------------------------------------------------------

## Save one panel as SVG (editable, used for the final figure assembly) and PNG
## (quick look).  `name` is the file stem, e.g. "Fig1B_Evenness".
save_panel <- function(plot, name, width, height, dpi = 300) {
  svg_file <- file.path(paths$figures, paste0(name, ".svg"))
  png_file <- file.path(paths$figures, paste0(name, ".png"))
  ggsave(svg_file, plot, width = width, height = height, device = svglite::svglite, bg = "white")
  ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white")
  message("saved ", svg_file)
  invisible(svg_file)
}

## Open an SVG + PNG device pair for grid-based graphics (ComplexHeatmap, igraph).
## Usage:  dev <- open_panel("Fig2C_KO_heatmap", 8, 6); draw(ht); close_panel(dev)
open_panel <- function(name, width, height, dpi = 300) {
  files <- c(svg = file.path(paths$figures, paste0(name, ".svg")),
             png = file.path(paths$figures, paste0(name, ".png")))
  list(files = files, width = width, height = height, dpi = dpi)
}
draw_panel <- function(dev, expr) {
  expr <- substitute(expr)
  env <- parent.frame()
  svglite::svglite(dev$files[["svg"]], width = dev$width, height = dev$height, bg = "white")
  eval(expr, env); grDevices::dev.off()
  ## type = "cairo": the default Windows bitmap device cannot draw non-ASCII
  ## characters such as the en dash or Greek letters in gene-set labels
  png_args <- list(dev$files[["png"]], width = dev$width, height = dev$height, units = "in", res = dev$dpi, bg = "white")
  if (isTRUE(capabilities("cairo"))) png_args$type <- "cairo"
  do.call(grDevices::png, png_args)
  eval(expr, env); grDevices::dev.off()
  message("saved ", dev$files[["svg"]])
  invisible(dev$files)
}

## Write a results table next to the figures (CSV, UTF-8)
save_table <- function(x, name) {
  f <- file.path(paths$tables, paste0(name, ".csv"))
  write.csv(x, f, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(f)
}

## ---- 7. Session record --------------------------------------------------------

## Print R and key package versions into the log of every run
report_versions <- function(pkgs) {
  ip <- utils::installed.packages()[, "Version"]
  v <- vapply(pkgs, function(p) if (p %in% names(ip)) ip[[p]] else "not installed", character(1))
  message(R.version.string, " | ", paste(paste0(pkgs, " ", v), collapse = ", "))
}
