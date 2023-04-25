library(tidyverse)
library(anota2seq)
library(biomaRt)
################################PARAMETERS######################################
filter_samples = "17AC" # NULL | 17AC | 02136
filter_genes = "custom" # custom | allzeros | NULL
exclude_samples = NULL # NULL c("Batch_A_17AC_FAKi", "Batch_A_17AC_TGF")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/07_stimulated_CAFs/01_17AC_Conditions/")

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

BatchVec <- dplyr::filter(manip_info, Fraction=="F8") %>%
  dplyr::select(Batch) %>% 
  deframe()

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
                          generateSingleGenePlots = TRUE)

ads <- anota2seqResidOutlierTest(ads)

## Analysis
phenoLev <- unique(phenoVec)
myContrast <- matrix(nrow =length(phenoLev),ncol=length(phenoLev)-1)
rownames(myContrast) <- phenoLev
myContrast[,1] = c(1,0,0,-1,0)
myContrast[,2]= c(0,1,0,-1,0)
myContrast[,3]= c(0,0,1,-1,0)
myContrast[,4]= c(0,0,0,-1,1)


ads <- anota2seqAnalyze(ads, contrasts = myContrast,
                 correctionMethod = "BH", useProgBar = TRUE, fileStem = "ANOTA2SEQ",
                 analysis = c("translation", "buffering", "translated mRNA", "total mRNA"))


#### Check threshold passing in contrast 1 (POSSIBLE ERROR)
ads@translation@apvStatsRvm[[3]] %>%
  as_tibble(rownames = "EnsemblID") %>% summary()

ads@buffering@apvStatsRvm[[3]] %>%
  as_tibble(rownames = "EnsemblID") %>% summary()

ads@translatedmRNA@apvStatsRvm[[3]] %>%
  as_tibble(rownames = "EnsemblID") %>% summary()

ads@totalmRNA@apvStatsRvm[[3]] %>%
  as_tibble(rownames = "EnsemblID") %>% summary()


ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            selContrast = 1,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 0.25)

ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            selContrast = 2,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 1)
## NEAA vs NT
ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            selContrast = 3,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 0.25)

anota2seqGetOutput(object = ads, output="regModes",
                   selContrast = 3, analysis="translation",
                   getRVM = TRUE)

anota2seqGetOutput(object = ads, output="regModes",
                   selContrast = 3, analysis="buffering",
                   getRVM = TRUE)

anota2seqGetOutput(object = ads, output="regModes",
                   selContrast = 3, analysis="mRNA abundance",
                   getRVM = TRUE)

ads <- anota2seqSelSigGenes(Anota2seqDataSet = ads,
                            selContrast = 4,
                            minSlopeTranslation = -1,
                            maxSlopeTranslation = 2,
                            minSlopeBuffering = -2,
                            maxSlopeBuffering = 1,
                            maxPAdj = 0.25)
anota2seqGetOutput(object = ads, output="regModes",
                   selContrast = 3, analysis="mRNA abundance",
                   getRVM = TRUE)

for (i in seq(1,4)){
  name = paste0(rownames(myContrast)[myContrast[,i]==1], " vs. NT")
  par(mfrow = c(1, 2))
  anota2seqPlotPvalues(ads, selContrast = i, useRVM =TRUE, plotToFile = FALSE, contrastName = name)

  ads <- anota2seqRegModes(ads, c(TRUE, TRUE))
  anota2seqPlotFC(ads, selContrast = i, plotToFile = FALSE, contrastName = name)
}

