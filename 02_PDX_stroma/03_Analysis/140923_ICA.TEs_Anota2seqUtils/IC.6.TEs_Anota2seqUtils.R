#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
library(ggVennDiagram)
library(msigdbr)
library(pheatmap)
library(Hmisc)
library(ggrepel)
source("../100122_ICABoot/functions.R")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/140923_ICA.TEs_Anota2seqUtils/")

# PARAMETERS
#-------------------------------------------------------------------------------
component_reorientation = TRUE
reorient_TEs <- c(1, 1, -1, 1, 1, -1)
#-------------------------------------------------------------------------------

# load annotation with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "mmusculus_gene_ensembl",
                        version = 75)

#listAttributes(ensembl75, page="feature_page")
annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id',
                                        'entrezgene',
                                        'mgi_id',
                                        'chromosome_name'), mart = ensembl75)

gene_to_ensembl = deframe(annot_ensembl75[c( "external_gene_id", "ensembl_gene_id")])
ensembl_to_gene = setNames(names(gene_to_ensembl), gene_to_ensembl)

# Data loading
## ICs
ICA_TEs <- read_rds("../100122_ICABoot/02_Output/ICA_TEs.RDS")

if (component_reorientation == TRUE){
  ICA_TEs[["A"]] <- t(t(ICA_TEs[["A"]]) * reorient_TEs)
  ICA_TEs[["S"]] <- t(t(ICA_TEs[["S"]]) * reorient_TEs)
}

## Norm_data
TEs_m <- read_tsv("../081221_TranslationEfficacy/02_Output/TEs.tsv")

## Metadata 
annotations <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>%
  dplyr::select(sample, PAMG, ICA1, ICA3) %>%
  dplyr::rename(ISRact = ICA3)

# Anota2seqUtils of the TE of IC.6 high samples
A_TEs <- as_tibble(ICA_TEs$A, rownames = "sample")

## selection of high and low samples
hard_threshold <- if_else(A_TEs$IC.6 > 0.1, "high", "low")
soft_threshold <- if_else(A_TEs$IC.6 > 0, "high", "low")

threshold = hard_threshold
A_TEs$threshold = threshold

ggplot(A_TEs, aes(IC.6)) +
  geom_density() +
  geom_rug(aes(color = threshold)) +
  theme_bw()

## data processing
TEs_subset <- dplyr::mutate(TEs_m, geneID = ensembl_to_gene[EnsemblID]) %>%
  distinct(geneID, .keep_all = TRUE) %>%
  drop_na(geneID) %>%
  dplyr::select(geneID, A_TEs$sample[threshold == "high"])

TEs_subset_mean <- dplyr::mutate(TEs_subset, IC.6high = rowMeans(TEs_subset[-1])) %>%
  dplyr::select(geneID, IC.6high) %>% 
  dplyr::arrange(abs(IC.6high))

TEs_5perc <- tail(TEs_subset_mean, n = 1000)
TEs_5perc_up <- dplyr::filter(TEs_5perc, IC.6high > 0)
TEs_5perc_down <- dplyr::filter(TEs_5perc, IC.6high < 0)
  
TEs_5perc_l <- list(translationUp = TEs_5perc_up$geneID,
                   translationDown = TEs_5perc_down$geneID)
## anota2seqUtils
library(anota2seqUtils)

annot <- retrieveFormatData(source = "load", species = "mouse")

### length
len_perc <- lengthAnalysis(geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
                      customBg = TEs_subset$geneID, # like so we need to give some background
                      regulation = c("translationUp", "translationDown"),
                      contrast = c(1,1), 
                      region = c('UTR5', 'CDS', 'UTR3'),
                      selection = 'longest', 
                      annot = annot,
                      plotType = 'boxplot',
                      pdfName = "results/extremeTEs/5percabs")


### Compare difference in nucleotide composition
content_perc <- contentAnalysis(geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
                           customBg = TEs_subset$geneID, # like so we need to give some background
                           regulation = c("translationUp", "translationDown"),
                           contrast = c(1,1), 
                           region = c('UTR5', 'CDS', 'UTR3'),
                           selection = 'longest', 
                           annot = annot, 
                           comparisons = list(c(0,1),c(0,2),c(1,2)),
                           contentIn = c("G", "C", "A", "T"),
                           plotType = 'boxplot',
                           pdfName = "results/extremeTEs/5percabs")

### Compare presence of uORFs
uorf_strong_perc <- uorf_analysis(geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
                             customBg = TEs_subset$geneID, # like so we need to give some background
                             regulation = c("translationUp", "translationDown"),
                             contrast = c(1,1),
                             onlyUTR5 = F,
                             startCodon = "ATG",
                             KozakContext = "strong",
                             selection = "longest",
                             annot = annot,
                             unitOut = "number",
                             pdfName = "results/extremeTEs/5percabs")

# ###  Run de novo sequence analysis
# deNovo <- motifAnalysis(annot = annot,
#                         memePath = "/home/jacobo/meme/bin/",
#                         geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
#                         customBg = TEs_subset$geneID, # like so we need to give some background
#                         regulation = c("translationUp", "translationDown"),
#                         contrast = c(1,1), 
#                         region = c('UTR5', 'CDS', 'UTR3'),
#                         subregion = NULL,
#                         subregionSel=NULL)
# 
# ###  Quantify presence of motifs
# motifs <- contentMotifs(annot = annot,
#                         geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
#                         customBg = TEs_subset$geneID, # like so we need to give some background
#                         regulation = c("translationUp", "translationDown"),
#                         contrast = c(1,1),
#                         motifsIn = deNovo[[1]]$motifsOut,
#                         region = "UTR3",
#                         dist = 1,
#                         selection = "longest",
#                         annot = annot,
#                         comparisons = list(c(1,2)),
#                         unitOut = "number",
#                         pdfName = "results/extremeTEs/5percabs"
# )

### Folding Energies
feOut_perc <- foldingEnergyAnalysis(geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
                               customBg = TEs_subset$geneID, # like so we need to give some background
                               species = 'mouse',
                               regulation = c("translationUp", "translationDown"),
                               contrast = c(1,1), 
                               region = c('UTR5', 'CDS', 'UTR3'),
                               selection = 'longuest', 
                               annot = annot,
                               plotType = 'ecdf',
                               pdfName = "results/extremeTEs/5percabs",
                               residFE = TRUE,
                               plotOut = TRUE)
# ### Codon Analysis (FAIL)
# codonOut_perc <- codonUsage(geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
#                        customBg = TEs_subset$geneID, # like so we need to give some background
#                        analysis = "codon",
#                        annot = annot,
#                        regulation = c("translationUp", "translationDown"),
#                        contrast = c(1,1),
#                        type = "sequence",
#                        comparisons = list(c(1,2)),
#                        annotType = 'ccds',
#                        sourceCod = "load",
#                        selection = "longest",
#                        pAdj=0.01,
#                        plotHeatmap = TRUE,
#                        pdfName = "results/extremeTEs/5percabs",
#                        species = "mouse")

# ### selCodonOut_perc <- codonCalc(codonIn = codonOut[['codonAll']], analysis = "codon",
#                          featsel = codonOut[[1]], 
#                          unit = "freq", 
#                          regulation = c("translationUp", "translationDown"),
#                          contrast = c(1,1),
#                          geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
#                          customBg = TEs_subset$geneID, # like so we need to give some background
#                          pdfName = 'Up',
#                          plotOut = T,
#                          plotType = 'ecdf')



### No human signatures, thus not straightforward testing

## Feature Integration
features_perc <- c(len_perc,
              content_perc,
              uorf_strong_perc,
              #motifs_perc,
              feOut_perc #,
              #selCodonOut, 
              #sign
              )

featureIntegration(geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
                   customBg = TEs_subset$geneID, # like so we need to give some background
                   effectMeasure = deframe(TEs_5perc),
                   contrastSel = 1,
                   features = features_perc,
                   pdfName = "results/extremeTEs/5percabs",
                   regOnly = T,
                   allFeat = F,
                   regulationGen = "translation",
                   analysis_type = "lm")

################################################################################
# Anota2seqUtils of the Gene weights
S_TEs <- as_tibble(ICA_TEs$S, rownames = "EnsemblID") %>%
  dplyr::mutate(geneID = ensembl_to_gene[EnsemblID]) %>%
  drop_na(geneID) %>% 
  distinct(geneID, .keep_all = TRUE) %>%
  dplyr::select(geneID, matches("IC."))

## selection of high and low genes
hard_threshold <- if_else(abs(S_TEs$IC.6) > 1, if_else(S_TEs$IC.6 > 0, "high", "low"), NA)
soft_threshold <- if_else(S_TEs$IC.6 > 0, "high", "low")

threshold = hard_threshold
S_TEs$threshold = threshold

ggplot(S_TEs, aes(IC.6)) +
  geom_density() +
  geom_rug(aes(color = threshold)) +
  theme_bw()

## data processing
S_TEs_extup <- dplyr::filter(S_TEs, threshold == "high")
S_TEs_extdown <- dplyr::filter(S_TEs, threshold == "low")

S_TEs_ext_l <- list(translationUp = S_TEs_extup$geneID,
                    translationDown = S_TEs_extdown$geneID)

### length
len_ext <- lengthAnalysis(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                      customBg = S_TEs$geneID, # like so we need to give some background
                      regulation = c("translationUp", "translationDown"),
                      contrast = c(1,1), 
                      region = c('UTR5', 'CDS', 'UTR3'),
                      selection = 'longest', 
                      annot = annot,
                      plotType = 'boxplot',
                      pdfName = "results/extremeIC.6GW/abs1")


### Compare difference in nucleotide composition
content_ext <- contentAnalysis(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                           customBg = S_TEs$geneID, # like so we need to give some background
                           regulation = c("translationUp", "translationDown"),
                           contrast = c(1,1), 
                           region = c('UTR5', 'CDS', 'UTR3'),
                           selection = 'longest', 
                           annot = annot, 
                           comparisons = list(c(0,1),c(0,2),c(1,2)),
                           contentIn = c("G", "C", "A", "T"),
                           plotType = 'boxplot',
                           pdfName = "results/extremeIC.6GW/abs1")

### Compare presence of uORFs
uorf_strong_ext <- uorf_analysis(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                             customBg = S_TEs$geneID, # like so we need to give some background
                             regulation = c("translationUp", "translationDown"),
                             contrast = c(1,1),
                             onlyUTR5 = F,
                             startCodon = "ATG",
                             KozakContext = "strong",
                             selection = "longest",
                             annot = annot,
                             unitOut = "number",
                             pdfName = "results/extremeIC.6GW/abs1")

###  Run de novo sequence analysis #BUG?
# deNovo <- motifAnalysis( geneList = S_TEs_ext_l, 
#                          customBg = S_TEs$geneID, 
#                          annot = annot,
#                          memePath = "/home/jacobo/meme/bin/",
#                          regulation = c("translationUp", "translationDown"),
#                          contrast = c(1,1),
#                          region = c('UTR5', 'CDS', 'UTR3'),
#                          subregion = NULL,
#                          subregionSel=NULL)
#
# ###  Quantify presence of motifs #BUG?
# motifs <- contentMotifs(annot = annot,
#                         geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
#                         customBg = S_TEs$geneID, # like so we need to give some background
#                         regulation = c("translationUp", "translationDown"),
#                         contrast = c(1,1),
#                         motifsIn = deNovo[[1]]$motifsOut,
#                         region = "UTR3",
#                         dist = 1,
#                         selection = "longest",
#                         annot = annot,
#                         comparisons = list(c(1,2)),
#                         unitOut = "number",
#                         pdfName = "results/extremeIC.6GW/abs1"
# )

### Folding Energies
feOut_ext <- foldingEnergyAnalysis(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                                    customBg = S_TEs$geneID, # like so we need to give some background
                                    species = 'mouse',
                                    regulation = c("translationUp", "translationDown"),
                                    contrast = c(1,1), 
                                    region = c('UTR5', 'CDS', 'UTR3'),
                                    selection = 'longuest', 
                                    annot = annot,
                                    plotType = 'ecdf',
                                    pdfName = "results/extremeIC.6GW/abs1",
                                    residFE = TRUE,
                                    plotOut = TRUE)

# ### Codon Analysis (FAIL)
# codonOut_perc <- codonUsage(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
#                             customBg = S_TEs$geneID, # like so we need to give some background
#                             analysis = "codon",
#                             annot = annot,
#                             regulation = c("translationUp", "translationDown"),
#                             contrast = c(1,1),
#                             type = "sequence",
#                             comparisons = list(c(1,2)),
#                             annotType = 'ccds',
#                             sourceCod = "load",
#                             selection = "longest",
#                             pAdj=0.01,
#                             plotHeatmap = TRUE,
#                             pdfName = "results/extremeIC.6GW/abs1",
#                             species = "mouse")

# ### selCodonOut_perc <- codonCalc(codonIn = codonOut[['codonAll']], analysis = "codon",
#                          featsel = codonOut[[1]], 
#                          unit = "freq", 
#                          regulation = c("translationUp", "translationDown"),
#                          contrast = c(1,1),
#                          geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
#                          customBg = TEs_subset$geneID, # like so we need to give some background
#                          pdfName = 'Up',
#                          plotOut = T,
#                          plotType = 'ecdf')



### No human signatures, thus not straightforward testing

## Feature Integration
features_ext <- c(len_ext,
                   content_ext,
                   uorf_strong_ext,
                   #motifs_perc,
                   feOut_ext #,
                   #selCodonOut, 
                   #sign
)

featureIntegration(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                   customBg = S_TEs$geneID, # like so we need to give some background
                   effectMeasure = dplyr::select(S_TEs,geneID, IC.6) %>% deframe(),
                   contrastSel = 1,
                   features = features_ext,
                   pdfName = "results/extremeIC.6GW/abs1",
                   regOnly = T,
                   allFeat = F,
                   regulationGen = "translation",
                   analysis_type = "lm")
