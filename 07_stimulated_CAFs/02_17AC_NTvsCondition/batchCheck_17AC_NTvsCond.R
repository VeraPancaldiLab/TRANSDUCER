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
phenoVec <- dplyr::filter(manip_info, Fraction=="F8") %>% # just to select one of the fractions
  dplyr::select(Condition) %>% 
  deframe()

lfc <- inner_join(dataP, dataT) %>%
  column_to_rownames("Geneid") %>% 
  DGEList() %>% 
  calcNormFactors(method = "TMM") %>%
  cpm() %>% 
  as_tibble(rownames = "Geneid") %>% 
  mutate(lfc_Input = log2(Batch_A_17AC_NEAA_Input/Batch_A_17AC_NT_Input ),
         lfc_F8 = log2(Batch_A_17AC_NEAA_F8/Batch_A_17AC_NT_F8 )) 

# add gene names
## get gene ID with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)

#listAttributes(ensembl75, page="feature_page")
annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)

translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

lfc <- dplyr::mutate(lfc, identifier = translate[Geneid])

# lfc scatter plot
stimuli_markers <- read_tsv(file = "../00_Data/stimuli_markers.tsv")
highlight_list <- dplyr::filter(stimuli_markers, stimuli %in% phenoVec) %>% deframe()

tibble(lfc) %>%
  mutate(highlight = if_else(identifier %in% highlight_list, identifier, "Other") %>% fct(levels = c("Other", highlight_list))) %>% 
  dplyr::arrange(highlight) %>%
  ggplot() +
  aes(x = lfc_Input, y = lfc_F8, color = highlight) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  theme_classic() + 
  labs(title = filter_samples)


