library(tidyverse)
library(anota2seq)
library(biomaRt)
################################PARAMETERS######################################
filter_samples = "A_17AC_(NT|NEAA)" # (A|B|C)_17AC_(NT|NEAA) | 17AC_(NT|IL1) | (A|B|D)_17AC_(NT|TGF) | 17AC_(NT|FAKi)
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

# Foldchange calculation
dataP <-  dplyr::select(counts, Geneid, deframe(manip_info[manip_info$Fraction=="F8","sample_name"]))
dataT <-  dplyr::select(counts, Geneid, deframe(manip_info[manip_info$Fraction=="Input","sample_name"]))

lfc <- inner_join(dataP, dataT) %>% 
  mutate(lfc_Input = Batch_A_17AC_NEAA_Input/Batch_A_17AC_NT_Input -1,
         lfc_F8 = Batch_A_17AC_NEAA_F8/Batch_A_17AC_NT_F8 -1) %>%
  dplyr::select(Geneid,lfc_Input,lfc_F8)

