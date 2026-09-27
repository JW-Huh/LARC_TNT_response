## =============================================================================
##  R/_helpers_enrichment_map.R  --  text rules used by 12_Fig4A_Host_enrichment_map.R
##  ---------------------------------------------------------------------------
##  Three keyword-based helpers that turn MSigDB gene-set identifiers into
##  readable labels, assign each gene set to a broad "biological axis", and
##  name a network community from the gene sets it contains.  They are only
##  needed when the enrichment map is recomputed (refit_gsea <- TRUE in script
##  12); in replay mode the frozen node table already carries these labels.
##  All rules are plain regular expressions - edit them to change the wording.
## =============================================================================

suppressPackageStartupMessages(library(stringr))

## "HALLMARK_G2M_CHECKPOINT" -> "G2/M checkpoint"
format_pathway_label <- function(x) {
  x %>%
    str_remove("^(HALLMARK|REACTOME|GOBP|GOCC|GOMF|WP|KEGG)_") %>%
    str_replace_all("_", " ") %>% str_to_lower() %>% str_to_sentence() %>%
    str_replace_all("\\bNfkb\\b", "NF-kB") %>% str_replace_all("\\bTnf\\b", "TNF") %>%
    str_replace_all("\\bIfn\\b", "IFN") %>% str_replace_all("\\bIl6\\b", "IL-6") %>%
    str_replace_all("\\bJak\\b", "JAK") %>% str_replace_all("\\bStat\\b", "STAT") %>%
    str_replace_all("\\bTlr\\b", "TLR") %>% str_replace_all("\\bNod\\b", "NOD") %>%
    str_replace_all("\\bAhr\\b", "AhR") %>% str_replace_all("\\bNad\\b", "NAD") %>%
    str_replace_all("\\bEmt\\b", "EMT") %>% str_replace_all("\\bRos\\b", "ROS") %>%
    str_replace_all("\\bG2m\\b", "G2/M") %>% str_replace_all("\\bMyc\\b", "MYC") %>%
    str_squish()
}

## Broad biological axis of a gene set from its label (first matching rule wins)
assign_biological_axis <- function(label) {
  u <- str_to_upper(label)
  dplyr::case_when(
    str_detect(u, "COLORECTAL|ADENOMA|\\bCRC\\b|WNT|BETA.?CATENIN|\\bMYC(?:N|L)?\\b|ONCOGEN") ~ "Adenoma/CRC epithelial program",
    str_detect(u, "G2.?M|E2F|CELL CYCLE|MITOTIC|DNA REPLICATION|CHROMOSOME SEGREGATION|CENTROSOME|SPINDLE|S PHASE|M PHASE") ~ "Cell-cycle checkpoint and proliferation",
    str_detect(u, "EPITHELIAL.?MESENCHYMAL|\\bEMT\\b|ECM|COLLAGEN|INTEGRIN|EXTRACELLULAR MATRIX|MATRIX|MATRISOME|FIBROSIS|FIBROTIC|MIGRATION|INVASION|ADHESION") ~ "Epithelial–mesenchymal and matrix remodeling",
    str_detect(u, "PROSTAGLANDIN|PROSTANOID|\\bPGE.?2\\b|\\bPTGS1\\b|\\bPTGS2\\b") ~ "Prostaglandin signaling and response",
    str_detect(u, "OXIDATIVE|\\bROS\\b|HYPOXIA|REACTIVE OXYGEN|NRF2|STRESS RESPONSE") ~ "Hypoxia and oxidative-stress response",
    str_detect(u, "T CELL|LYMPH|LEUKOCYTE|ANTIGEN|INTERFERON|\\bIFN\\b|IMMUNE|IMMUN|CYTOKINE|NF.?KB|\\bTNF[A-Z0-9]*\\b|IL.?6|JAK|\\bSTAT[0-9AB]*\\b|\\bTLR[0-9]*\\b|\\bNOD[0-9]*\\b|NOD.?LIKE") ~ "Immune receptor and cytokine signaling",
    str_detect(u, "MUCUS|GOBLET|MUCIN|BARRIER|TIGHT JUNCTION|APICAL|ANTIMICROBIAL|SECRETORY") ~ "Barrier and mucus program",
    str_detect(u, "\\bAHR\\b|INDOLE|TRYPTOPHAN|\\bNAD\\b|\\bNADH\\b|\\bNADP\\b|\\bNADPH\\b|NIACIN|NICOTIN|LIPID|PHOSPHOLIPID|GLYCOPROTEIN|GLYCAN|METABOL") ~ "Lipid/glycoprotein and metabolic remodeling",
    str_detect(u, "RIBOSOM|TRANSLATION|TRANSLATIONAL|RNA PROCESSING|MRNA|PROTEIN SYNTHESIS") ~ "Ribosome and translation initiation",
    str_detect(u, "GOLGI|VESICLE|ENDOSOME|MEMBRANE TRAFFICKING|COATED VESICLE|SECRETION") ~ "Golgi–endosomal vesicle trafficking",
    TRUE ~ "Other"
  )
}

## Name a community from the labels of its member gene sets.  The rules are
## ordered from specific (MYC, EMT, adenoma ...) to generic (fallback = the
## most frequent biological axis).
refine_cluster_label <- function(term_vec, axis_vec) {
  txt <- str_to_upper(paste(term_vec, collapse = " | "))
  terms <- str_to_upper(term_vec)
  n_hyp  <- sum(str_detect(terms, "HYPOXIA"))
  n_oxi  <- sum(str_detect(terms, "OXIDATIVE|\\bROS\\b|REACTIVE OXYGEN") & !str_detect(terms, "HYPOXIA"))
  n_pg   <- sum(str_detect(terms, "PROSTAGLANDIN|PROSTANOID|\\bPGE.?2\\b|\\bPTGS1\\b|\\bPTGS2\\b"))

  myc_up <- str_detect(txt, "\\bMYC\\b.*\\bUP\\b|\\bUP\\b.*\\bMYC\\b|MYC TARGET|HALLMARK MYC")
  myc_dn <- str_detect(txt, "\\bMYC\\b.*\\bDN\\b|\\bDN\\b.*\\bMYC\\b|MYC.*DOWN|DOWN.*MYC")
  if (myc_up && !myc_dn) return("MYC target\nprogram")
  if (myc_dn && !myc_up) return("MYC-down\nsignature")
  if (myc_up && myc_dn)  return("Mixed MYC\nCRC signatures")
  if (str_detect(txt, "EPITHELIAL.?MESENCHYMAL|\\bEMT\\b")) return("Epithelial–mesenchymal\ntransition")
  if (str_detect(txt, "ECM|COLLAGEN|EXTRACELLULAR MATRIX|MATRIX|MATRISOME|FIBROSIS|FIBROTIC")) return("ECM/collagen\nremodeling")
  if (str_detect(txt, "ADENOMA.*UP")) return("Colorectal adenoma\nsignature")
  if (str_detect(txt, "ADENOMA.*DN")) return("Colorectal adenoma-down\nsignature")
  if (str_detect(txt, "ADENOMA")) return("Colorectal adenoma\nsignature")
  if (str_detect(txt, "WNT|BETA.?CATENIN")) return("Wnt/beta-catenin\nCRC program")
  if (str_detect(txt, "\\bMYC\\b")) return("MYC-associated\nCRC signature")
  if (str_detect(txt, "G2.?M")) return("G2/M checkpoint")
  if (str_detect(txt, "E2F")) return("E2F target\nprogram")
  if (str_detect(txt, "CELL CYCLE|MITOTIC|DNA REPLICATION|PROLIFERATION|CENTROSOME|SPINDLE")) return("Cell-cycle\nprogression")

  fam <- c(hypoxia = n_hyp, oxidative = n_oxi, prostaglandin = n_pg)
  if (max(fam) > 0) {
    top <- names(fam)[fam == max(fam)]
    if (length(top) == 1) return(switch(top, hypoxia = "Hypoxia\nresponse", oxidative = "Oxidative-stress /\nROS response",
                                        prostaglandin = "Prostaglandin\nresponse"))
    if (setequal(top, c("hypoxia", "oxidative"))) return("Hypoxia / redox-\nrelated programs")
    if (setequal(top, c("hypoxia", "prostaglandin"))) return("Hypoxia / prostaglandin\nresponse")
    if (setequal(top, c("oxidative", "prostaglandin"))) return("Prostaglandin / redox\nresponse")
    return("Mixed stress-response\nprograms")
  }
  if (str_detect(txt, "\\bTLR[0-9]*\\b|\\bNOD[0-9]*\\b|NOD.?LIKE|PATTERN RECOGNITION")) return("TLR/NOD pattern-\nrecognition signaling")
  if (str_detect(txt, "\\bTNF[A-Z0-9]*\\b|NF.?KB|IL.?6|JAK|\\bSTAT[0-9AB]*\\b")) return("TNF–NFkB /\nIL-6–JAK–STAT signaling")
  if (str_detect(txt, "CHEMOKINE")) return("Chemokine\nsignaling")
  if (str_detect(txt, "CYTOKINE|INTERLEUKIN")) return("Cytokine\nsignaling")
  if (str_detect(txt, "T CELL|T-CELL|LYMPHOCYTE")) return("T-cell / lymphocyte\nactivation")
  if (str_detect(txt, "INTERFERON|\\bIFN\\b")) return("Interferon\nresponse")
  if (str_detect(txt, "IMMUNE|IMMUN|ANTIGEN")) return("Immune receptor\nsignaling")
  if (str_detect(txt, "MUCUS|GOBLET|MUCIN")) return("Mucus/goblet-cell\nprogram")
  if (str_detect(txt, "BARRIER|TIGHT JUNCTION|APICAL")) return("Barrier / apical\njunction program")
  if (str_detect(txt, "\\bNAD\\b|\\bNADH\\b|\\bNADP\\b|\\bNADPH\\b|NIACIN|NICOTIN")) return("NAD/niacin\nmetabolism")
  if (str_detect(txt, "\\bAHR\\b|INDOLE")) return("AhR/indole-ligand\nresponse")
  if (str_detect(txt, "LIPID|PHOSPHOLIPID|GLYCOPROTEIN|GLYCAN|METABOL")) return("Lipid/glycoprotein\nmetabolic remodeling")
  if (str_detect(txt, "RIBOSOM")) return("Ribosome-associated\ntranslation")
  if (str_detect(txt, "TRANSLATION|TRANSLATIONAL")) return("Translation\ninitiation")
  if (str_detect(txt, "GOLGI|VESICLE|ENDOSOME|TRAFFICKING")) return("Golgi–endosomal\nvesicle trafficking")
  if (str_detect(txt, "COLORECTAL|COLON|RECTAL|\\bCRC\\b")) return("CRC epithelial\nprogram")

  top_axis <- names(sort(table(axis_vec), decreasing = TRUE))[1]
  c("Adenoma/CRC epithelial program" = "CRC epithelial\nprogram",
    "Cell-cycle checkpoint and proliferation" = "Cell-cycle\nprogram",
    "Epithelial–mesenchymal and matrix remodeling" = "Epithelial–mesenchymal\ntransition",
    "Prostaglandin signaling and response" = "Prostaglandin\nresponse",
    "Hypoxia and oxidative-stress response" = "Hypoxia / oxidative-\nstress programs",
    "Immune receptor and cytokine signaling" = "Immune receptor\nsignaling",
    "Barrier and mucus program" = "Barrier / mucus\nprogram",
    "Lipid/glycoprotein and metabolic remodeling" = "Metabolic\nremodeling",
    "Ribosome and translation initiation" = "Ribosome-associated\ntranslation",
    "Golgi–endosomal vesicle trafficking" = "Golgi–endosomal\nvesicle trafficking")[top_axis] %||% top_axis
}
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
