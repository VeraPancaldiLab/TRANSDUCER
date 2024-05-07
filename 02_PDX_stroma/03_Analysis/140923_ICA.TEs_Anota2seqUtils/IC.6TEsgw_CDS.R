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
library(anota2seqUtils)

source("../100122_ICABoot/functions.R")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/140923_ICA.TEs_Anota2seqUtils/")

# PARAMETERS
#-------------------------------------------------------------------------------
component_reorientation = TRUE
reorient_TEs <- c(1, 1, -1, 1, 1, -1)
pdf_name <- "results/CDS/IC.6gw_"
# motif analysis
de_novo = FALSE
min_len = 5 #filter for min len
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


## Compare difference in nucleotide composition
content <- contentAnalysis(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                               customBg = S_TEs$geneID, # like so we need to give some background
                               geneListcolours = brewer.pal(n = 3, name = "RdBu")[c(1,3)],
                               regulation = c("translationUp", "translationDown"),
                               contrast = c(1,1), 
                               region = c('CDS'),
                               selection = 'longest', 
                               annot = annot, 
                               comparisons = list(c(0,1),c(0,2),c(1,2)),
                               contentIn = c("G", "C", "A", "T"),
                               plotType = 'boxplot',
                               pdfName = pdf_name)

## Motif Quantification
### Run de novo sequence analysis #BUG?
if (de_novo == TRUE){
motif_list <- motifAnalysis( geneList = S_TEs_ext_l,
                         customBg = S_TEs$geneID,
                         annot = annot,
                         memePath = "/home/jacobo/meme/bin/",
                         regulation = c("translationUp", "translationDown"),
                         contrast = c(1,1),
                         region = c('CDS'),
                         subregion = NULL,
                         subregionSel=NULL)
} else {
  motif_list <- c("DRACH", "HCARD")
  attract <- read_tsv("Input/ATtRACT_db.tsv") %>%
    filter(Organism == "Mus_musculus",
           Len >= min_len) %>% 
    dplyr::select(Gene_name, Motif) %>%
    deframe()
  motifs_in = motif_list # attract | motif_list
}


### Quantification
motifs_quant <- contentMotifs(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                            customBg = S_TEs$geneID, # like so we need to give some background
                            geneListcolours =  brewer.pal(n = 3, name = "RdBu")[c(1,3)],
                            regulation = c("translationUp", "translationDown"),
                            contrast = c(1,1),
                            motifsIn = motifs_in, 
                            region = c("CDS"),
                            dist = 1,
                            selection = "longest",
                            annot = annot,
                            comparisons = list(c(1,2)),
                            unitOut = "number",
                            pdfName = pdf_name)

## Codon Analysis # To be fixed
### caclulate usage per transcript
codonOut <- codonUsage(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                            customBg = S_TEs$geneID, # like so we need to give some background
                            geneListcolours =  brewer.pal(n = 3, name = "RdBu")[c(1,3)],
                            analysis = "codon",
                            annot = annot,
                            regulation = c("translationUp", "translationDown"),
                            contrast = c(1,1),
                            type = "sequence",
                            comparisons = list(c(1,2)),
                            annotType = 'ccds',
                            sourceCod = "load",
                            selection = "longest",
                            pAdj=0.01,
                            plotHeatmap = TRUE,
                            pdfName = pdf_name,
                            species = "mouse",)

### Statistical analysis
selCodonOut <- codonCalc(codonIn = codonOut[['codonAll']],
                          analysis = "codon",
                         featsel = codonOut[[1]],
                         geneListcolours =  brewer.pal(n = 3, name = "RdBu")[c(1,3)],
                         unit = "freq",
                         regulation = c("translationUp", "translationDown"),
                         contrast = c(1,1),
                         geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
                         customBg = TEs_subset$geneID, # like so we need to give some background
                         pdfName = pdf_name,
                         plotOut = T,
                         plotType = 'ecdf')

# Feature integration
features <- c(content,
              motifs_quant
              # selCodonOut, # To be fixed
)

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

