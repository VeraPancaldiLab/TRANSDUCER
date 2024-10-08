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
source("../100122_ICABoot/functions.R")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/140923_ICA.TEs_Anota2seqUtils/")

# Functions 
## Alternative plotting as fuction do not work
convert_to_numeric <- function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
}

plot_bar_enrich <- function(x, p.adjust_th, Title) {
  mutate(x, GeneRatioNum = sapply(GeneRatio, convert_to_numeric)) %>% 
    arrange(GeneRatioNum) %>% 
    mutate(Description = as_factor(Description)) %>% 
    dplyr::filter(p.adjust < p.adjust_th) %>%
    slice_tail(n=20) %>%
    ggplot(aes(x = GeneRatioNum, y = Description, fill = p.adjust)) +
    geom_bar(stat="identity") +
    scale_fill_gradient(low = "red", high = "blue", guide=guide_colourbar(reverse = TRUE), limits = c(10^-11,0.1), trans = "log", breaks = c( 0.05, 0.01, 0.001, 10^-5, 10^-10),oob=squish)+    theme_bw() +
    labs(title = Title) +
    xlab("Gene Ratio") +
    ylab("")
}

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
                      geneListcolours = c(brewer.pal(8, "Reds")[c(4,8)]),
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
                           geneListcolours = c(brewer.pal(8, "Reds")[c(4,8)]),
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
                             geneListcolours = c(brewer.pal(8, "Reds")[c(4,8)]),
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
motifs_perc <- contentMotifs(geneList = TEs_5perc_l, #instead of anota2seq object  you input geneList and effect_measures
                        customBg = TEs_subset$geneID, # like so we need to give some background
                        geneListcolours = c(brewer.pal(8, "Reds")[c(4,8)]),
                        regulation = c("translationUp", "translationDown"),
                        contrast = c(1,1),
                        motifsIn = "DRACH", #deNovo[[1]]$motifsOut
                        region = c("UTR5", "CDS","UTR3"),
                        dist = 1,
                        selection = "longest",
                        annot = annot,
                        comparisons = list(c(1,2)),
                        unitOut = "number",
                        pdfName = "results/extremeTEs/5percabs_"
)

### Folding Energies
feOut_perc <- foldingEnergyAnalysis(geneList = TEs_5perc_l, 
                               customBg = TEs_subset$geneID,
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



### Mice signatures given by Inci
load("mouseSignatures.RData")

sign_perc <- signCalc(geneList = TEs_5perc_l,
                 customBg = TEs_subset$geneID,
                 addSign = mouseSignatures,
                 annot = annot)

## Feature Integration
features_perc <- c(len_perc,
              content_perc,
              uorf_strong_perc,
              motifs_perc,
              feOut_perc,
              #selCodonOut, 
              sign_perc
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

## Enrichment
### Reactome
TEs_5perc_up_enrich_PA <- enrichPathway(gene= as.character(na.exclude(gene_to_entrez[TEs_5perc_l$translationUp])),
                                     universe = as.character(na.exclude(gene_to_entrez[TEs_subset$geneID])),
                                            pvalueCutoff = 0.05, readable=TRUE, organism = "mouse")

TEs_5perc_down_enrich_PA <- enrichPathway(gene= as.character(na.exclude(gene_to_entrez[TEs_5perc_l$translationDown])),
                                       universe = as.character(na.exclude(gene_to_entrez[TEs_subset$geneID])),
                                              pvalueCutoff = 0.05, readable=TRUE, organism = "mouse")
#### plots
TEs_5perc_up_enrich_PA@result %>%  
  plot_bar_enrich(0.05,"Reactome Upregulated in IC.6 (extreme mean residuals)")

TEs_5perc_down_enrich_PA@result %>%  
  plot_bar_enrich(0.05,"Reactome Downregulated in IC.6 (extreme mean residuals)")

### GO
ont = "MF"
TEs_5perc_up_enrich_GO <- enrichGO(gene = as.character(na.exclude(gene_to_entrez[TEs_5perc_l$translationUp])),
         universe = as.character(na.exclude(gene_to_entrez[TEs_subset$geneID])),
         ont      = ont,
         pvalueCutoff = 0.05, readable=TRUE, OrgDb = org.Mm.eg.db)

TEs_5perc_down_enrich_GO <- enrichGO(gene = as.character(na.exclude(gene_to_entrez[TEs_5perc_l$translationDown])),
                                   universe = as.character(na.exclude(gene_to_entrez[TEs_subset$geneID])),
                                   ont      = ont,
                                   pvalueCutoff = 0.05, readable=TRUE, OrgDb = org.Mm.eg.db)

#### plots
TEs_5perc_up_enrich_GO@result %>%  
  plot_bar_enrich(0.05,paste0("GO-", ont," Upregulated in IC.6 (extreme mean residuals)"))

TEs_5perc_down_enrich_GO@result %>%  
  plot_bar_enrich(0.05, paste0("GO-", ont," Downregulated in IC.6 (extreme mean residuals)"))

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

ggsave("results/Figures/ic6tes_distrib.svg",
       width = 6,
       height = 2)


## data processing
S_TEs_extup <- dplyr::filter(S_TEs, threshold == "high")
S_TEs_extdown <- dplyr::filter(S_TEs, threshold == "low")

S_TEs_ext_l <- list(translationUp = S_TEs_extup$geneID,
                    translationDown = S_TEs_extdown$geneID)

### length
len_ext <- lengthAnalysis(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                      customBg = S_TEs$geneID, # like so we need to give some background
                      geneListcolours = c(brewer.pal(8, "Reds")[c(4,8)]),
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
                           geneListcolours = c(brewer.pal(8, "Reds")[c(4,8)]),
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
                             geneListcolours = c(brewer.pal(8, "Reds")[c(4,8)]),
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
###  Quantify presence of motifs #BUG?
motifs_ext <- contentMotifs(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                        customBg = S_TEs$geneID, # like so we need to give some background
                        geneListcolours = c(brewer.pal(8, "Reds")[c(4,8)]),
                        regulation = c("translationUp", "translationDown"),
                        contrast = c(1,1),
                        motifsIn = "DRACH", #deNovo[[1]]$motifsOut
                        region = c("UTR5", "CDS","UTR3"),
                        dist = 1,
                        selection = "longest",
                        annot = annot,
                        comparisons = list(c(1,2)),
                        unitOut = "number",
                        pdfName = "results/extremeIC.6GW/abs1"
)

### Folding Energies
feOut_ext <- foldingEnergyAnalysis(geneList = S_TEs_ext_l, #instead of anota2seq object  you input geneList and effect_measures
                                    customBg = S_TEs$geneID, # like so we need to give some background
                                    species = 'mouse',
                                   geneListcolours = c(brewer.pal(8, "Reds")[c(4,8)]),
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



### Mice signatures given by Inci
sign_ext <- signCalc(geneList = S_TEs_ext_l,
                      customBg = S_TEs$geneID,
                      addSign = mouseSignatures,
                      annot = annot)



## Feature Integration
features_ext <- c(len_ext,
                   content_ext,
                   uorf_strong_ext,
                   motifs_ext,
                   feOut_ext ,
                   #selCodonOut, 
                   sign_ext
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

## Enrichment
### Reactome
TEs_5perc_up_enrich_PA <- enrichPathway(gene= as.character(na.exclude(gene_to_entrez[S_TEs_ext_l$translationUp])),
                                        universe = as.character(na.exclude(gene_to_entrez[S_TEs$geneID])),
                                        pvalueCutoff = 0.05, readable=TRUE, organism = "mouse")

TEs_5perc_down_enrich_PA <- enrichPathway(gene= as.character(na.exclude(gene_to_entrez[S_TEs_ext_l$translationDown])),
                                          universe = as.character(na.exclude(gene_to_entrez[S_TEs$geneID])),
                                          pvalueCutoff = 0.05, readable=TRUE, organism = "mouse")
#### plots
reactome_pos_plot <- TEs_5perc_up_enrich_PA@result %>%  
  plot_bar_enrich(0.05,"Reactome Upregulated in IC.6 (extreme gene weights)")

reactome_neg_plot <- TEs_5perc_down_enrich_PA@result %>%  
  plot_bar_enrich(0.05,"Reactome Downregulated in IC.6 (extreme gene weights)")

### GO
pos_GO_plots = list()
neg_GO_plots = list()
for (ont in c("MF", "BP", "CC")){
  TEs_5perc_up_enrich_GO <- enrichGO(gene = as.character(na.exclude(gene_to_entrez[S_TEs_ext_l$translationUp])),
                                     universe = as.character(na.exclude(gene_to_entrez[S_TEs$geneID])),
                                     ont      = ont,
                                     pvalueCutoff = 0.05, readable=TRUE, OrgDb = org.Mm.eg.db)
  
  TEs_5perc_down_enrich_GO <- enrichGO(gene = as.character(na.exclude(gene_to_entrez[S_TEs_ext_l$translationDown])),
                                       universe = as.character(na.exclude(gene_to_entrez[S_TEs$geneID])),
                                       ont      = ont,
                                       pvalueCutoff = 0.05, readable=TRUE, OrgDb = org.Mm.eg.db)
  
  #### plots
  pos_GO_plots[[ont]] <- TEs_5perc_up_enrich_GO@result %>%  
    plot_bar_enrich(0.05,paste0("GO-", ont," Upregulated in IC.6 (extreme gene weights)"))
  
  neg_GO_plots[[ont]] <- TEs_5perc_down_enrich_GO@result %>%  
    plot_bar_enrich(0.05, paste0("GO-", ont," Downregulated in IC.6 (extreme gene weights)"))
}

# create joint positive plot for export
(reactome_pos_plot) / (pos_GO_plots[["MF"]]) / (pos_GO_plots[["BP"]]) / (pos_GO_plots[["CC"]]) + plot_layout(heights = c(1+5,
                                                                                                                         1+19,
                                                                                                                                1,
                                                                                                                                1), guides = 'collect')
ggsave("results/Figures/enrichments_pos.svg",
       width = 8,
       height = 7)

# create joint negative plot for export
(reactome_neg_plot) / (neg_GO_plots[["MF"]]) /(neg_GO_plots[["BP"]]) /(neg_GO_plots[["CC"]]) + plot_layout(heights = c(1+20,
                                                                                                                       1+20,
                                                                                                                       1+20,
                                                                                                                       1+20), guides = 'collect')
ggsave("results/Figures/enrichments_neg.svg",
       width = 8,
       height = 15)

