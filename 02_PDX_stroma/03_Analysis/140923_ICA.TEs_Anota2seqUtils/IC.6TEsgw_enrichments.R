#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
library(ggVennDiagram)
library(msigdbr)
library(pheatmap)
library(Hmisc)
library(ggrepel)
library(ReactomePA)
library(clusterProfiler)
library(RColorBrewer)
library(anota2seqUtils) # devtools::install('/home/jacobo/Anota2seq_handson/anota2sequtils-main')

source("../100122_ICABoot/functions.R")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/140923_ICA.TEs_Anota2seqUtils/")

# PARAMETERS
#-------------------------------------------------------------------------------
component_reorientation = TRUE
reorient_TEs <- c(1, 1, -1, 1, 1, -1)
pdf_name <- "results/enrichments/IC.6gw_"
de_novo = FALSE
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
gene_to_entrez = deframe(annot_ensembl75[c( "external_gene_id", "entrezgene")])
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

## Annotation
annot <- retrieveFormatData(source = "load", species = "mouse")

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

## Mice signatures given by Inci
load("mouseSignatures.RData")
signEnrch <- signCalc(geneList = S_TEs_ext_l,
                      customBg = S_TEs$geneID,
                      addSign = mouseSignatures,
                      annot = annot)

# Feature integration
features <- c(signEnrch)

featureIntegration(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                   customBg = S_TEs$geneID, # like so we need to give some background
                   effectMeasure = dplyr::select(S_TEs,geneID, IC.6) %>% deframe(),
                   contrastSel = 1,
                   features = features,
                   pdfName = pdf_name,
                   regOnly = T,
                   allFeat = F,
                   regulationGen = "translation",
                   analysis_type = "lm",
                   geneListcolours = brewer.pal(n = 3, name = "RdBu")[c(1,3)])
