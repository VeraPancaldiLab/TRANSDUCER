library(tidyverse)
library(anota2seq)
library(biomaRt)
################################PARAMETERS######################################
filter_samples = "17AC_(NT|TGF)" # (A|B|C)_17AC_(NT|NEAA) | 17AC_(NT|IL1) | (A|B|D)_17AC_(NT|TGF) | 17AC_(NT|FAKi)
filter_genes = "allzeros" # custom | allzeros | NULL
exclude_samples = NULL # NULL c("Batch_A_17AC_FAKi", "Batch_A_17AC_TGF")
correct_batch = TRUE
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/02_17AC_NTvsCondition/")

# Data loading
counts <- read_tsv("../00_Data/rawcounts.tsv")

## filter genes
if (filter_genes == "custom"){
  anotafilter = FALSE
  counts <- column_to_rownames(counts, "Geneid") %>%
    .[!(apply(., 1, function(x) {
      sum(x == 0) > ncol(counts)/2
    })), ] %>% as_tibble(rownames = "Geneid")
} else if (filter_genes == "allzeros"){
  anotafilter = TRUE
} else if (filter_genes == NULL){
  anotafilter = FALSE
}


## filter of samples
if (!is.null(filter_samples)){
  counts <- dplyr::select(counts, Geneid, names(counts)[str_detect(names(counts), pattern = filter_samples)])
}

if (!is.null(exclude_samples)){
  counts <- dplyr::select(counts, Geneid, !matches(exclude_samples))
}

## Metadata
manip_info <- read_tsv("../00_Data/Data_RNA_sample_Jacobo.tsv") %>% 
  dplyr::filter(sample_name %in% names(counts)) %>%
  dplyr::mutate(sample_name = fct(sample_name, levels=names(counts)[-1]),
                Condition = fct(Condition, levels= c("NT", unique(Condition)[ !unique(Condition) == 'NT']))) %>% # to reinforce NT as the first term always
  dplyr::arrange(sample_name)

all(names(counts)[-1] == manip_info$sample_name) %>% stopifnot()

# Anota2seq
dataP <-  dplyr::select(counts, Geneid, deframe(manip_info[manip_info$Fraction=="F8","sample_name"]))
dataT <-  dplyr::select(counts, Geneid, deframe(manip_info[manip_info$Fraction=="Input","sample_name"]))
phenoVec <- dplyr::filter(manip_info, Fraction=="F8") %>% # just to select one of the fractions
  dplyr::select(Condition) %>% 
  deframe()

if (correct_batch == TRUE){
  BatchVec <- dplyr::filter(manip_info, Fraction=="F8") %>%
  dplyr::select(Batch) %>% 
  deframe()
} else if (correct_batch == FALSE){
  BatchVec <- NULL
}

ads <- anota2seqDataSetFromMatrix(
  dataP = column_to_rownames(dataP, "Geneid"),
  dataT = column_to_rownames(dataT, "Geneid"),
  phenoVec = phenoVec,
  batchVec = BatchVec,
  dataType = "RNAseq",
  filterZeroGenes = anotafilter, # determined in filter_genes
  normalize = TRUE,
  transformation = "TMM-log2",
  varCutOff = NULL)

## QC
# ads <- anota2seqPerformQC(ads,
#                           generateSingleGenePlots = TRUE, 
#                           fileStem = "02_Output/")
# 
# ads <- anota2seqResidOutlierTest(ads)

## Analysis
### set contrast to reinforce correct order to NEAA vs NT
phenoLev <- unique(phenoVec)
myContrast <- matrix(nrow =length(phenoLev),ncol=length(phenoLev)-1)
rownames(myContrast) <- phenoLev
myContrast[,1] = c(-1, 1)

ads <- anota2seqAnalyze(ads,
                        correctionMethod = "BH", useProgBar = TRUE, fileStem = "02_Output/", contrasts = myContrast,
                        analysis = c("translation", "buffering", "translated mRNA", "total mRNA"))

ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            selContrast = 1,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 0.10)
## Result Plots
par(mfrow = c(1, 2))
anota2seqPlotPvalues(ads, selContrast = 1, useRVM = TRUE, plotToFile = FALSE, contrastName = paste(unique(phenoVec), collapse = " vs. "))

ads <- anota2seqRegModes(ads, c(TRUE, TRUE))
anota2seqPlotFC(ads, selContrast = 1, plotToFile = FALSE,  contrastName = paste(unique(phenoVec), collapse = " vs. "))

## Export results
ads_results <- anota2seqGetOutput(ads, output="singleDf", selContrast=1)


### get gene ID with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)

#listAttributes(ensembl75, page="feature_page")
annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)

translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

ads_results <- dplyr::mutate(ads_results, identifier = translate[identifier])

filename <- paste0("02_Output/result_tables/", paste(unique(phenoVec), collapse = "vs"), "_CorrectBatch", correct_batch, ".tsv")
  
write_tsv(ads_results, filename)

## Plot of specific factors/Markers
stimuli_markers <- read_tsv(file = "../00_Data/stimuli_markers.tsv")
highlight_list <- dplyr::filter(stimuli_markers, stimuli %in% phenoVec) %>% deframe()

marker_transplot <- tibble(ads_results) %>%
  mutate(highlight = if_else(identifier %in% highlight_list, identifier, "Other") %>% fct(levels = c("Other", highlight_list))) %>% 
  dplyr::arrange(highlight) %>%
  ggplot() +
  aes(x = totalmRNA.apvEff, y = translatedmRNA.apvEff, color = highlight) +
  scale_color_manual(values = c(Other = "grey", setNames(brewer.pal(n = length(highlight_list), name = "Set1"), c(highlight_list)))) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", colour = "red")+
  geom_hline(yintercept = 0,linetype = "dotted", colour = "red")+
  geom_vline(xintercept = 0, linetype = "dotted", colour = "red")+
  theme_classic() + 
  labs(title = filename)

ggsave(paste0("../FIGURES/marker_transplot_",paste(unique(phenoVec), collapse = "vs"), "_CorrectBatch", correct_batch, ".svg"),
       marker_transplot, width = 5,height = 4)


