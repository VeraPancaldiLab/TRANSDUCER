library(tidyverse)
library(anota2seq)
library(biomaRt)
################################PARAMETERS######################################
filter_samples = "(A|B|C)_17AC_(NT|NEAA)" # NULL | 17AC | 02136
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
  counts <- dplyr::select(counts, Geneid, names(counts)[!names(counts) %in% exclude_samples])
}

## Metadata
manip_info <- read_tsv("../00_Data/Data_RNA_sample_Jacobo.tsv") %>% 
  dplyr::filter(sample_name %in% names(counts)) %>%
  dplyr::mutate(sample_name = fct(sample_name, levels=names(counts)[-1])) %>%
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
ads <- anota2seqPerformQC(ads,
                          generateSingleGenePlots = TRUE, 
                          fileStem = "02_Output/")

ads <- anota2seqResidOutlierTest(ads)

## Analysis
ads <- anota2seqAnalyze(ads,
                        correctionMethod = "BH", useProgBar = TRUE, fileStem = "02_Output/",
                        analysis = c("translation", "buffering", "translated mRNA", "total mRNA"))

ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            selContrast = 1,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 0.25)
## Result Plots
par(mfrow = c(1, 2))
anota2seqPlotPvalues(ads, selContrast = 1, useRVM = TRUE, plotToFile = FALSE, contrastName = "provisional_name")

ads <- anota2seqRegModes(ads, c(TRUE, TRUE))
anota2seqPlotFC(ads, selContrast = 1, plotToFile = FALSE, contrastName = "WT vs. SRCinh")

