library(tidyverse)
library(anota2seq)
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
