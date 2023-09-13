library(tidyverse)
library(anota2seq)
library(anota2seqUtils)
library(biomaRt)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/081221_TranslationEfficacy/")

# Data loading
## Sample info and selection of data
include_med = FALSE

sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv")

top_samples <- dplyr::arrange(sample_info, ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low", "high"))

sample_info <- dplyr::select(top_samples, sample, ISRact) %>%
  right_join(sample_info, by= "sample") %>% 
  dplyr::mutate(ISRact = if_else(is.na(ISRact), "medium", ISRact))

if (include_med == FALSE){
  sample_info <- column_to_rownames(top_samples, "sample")
} else {
  sample_info <- column_to_rownames(sample_info, "sample")
}

# load and filter expression data
rowTumor_cyt <- read_tsv("../../00_Data/Processed_data/rawTumor_Cyt.tsv") %>% 
  dplyr::select(EnsemblID, all_of(rownames(sample_info)))
rowTumor_pol <- read_tsv("../../00_Data/Processed_data/rawTumor_Pol.tsv") %>%
  dplyr::select(EnsemblID, all_of(rownames(sample_info)))

# Translate to geneID
## get gene ID with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)

translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

rowTumor_cyt_gn <- dplyr::mutate(rowTumor_cyt, geneID = translate[EnsemblID]) %>%
  dplyr::select(geneID, matches("PDAC")) %>% drop_na(geneID) %>%
  distinct(geneID, .keep_all = T)

rowTumor_pol_gn <- dplyr::mutate(rowTumor_pol, geneID = translate[EnsemblID]) %>% 
  dplyr::select(geneID, matches("PDAC")) %>% drop_na(geneID) %>%
  distinct(geneID, .keep_all = T)

# Prepare data for anota2seq
## assert order
stopifnot(all(names(rowTumor_cyt_gn) == names(rowTumor_pol_gn)))
stopifnot(all(names(rowTumor_cyt_gn)[-1] == rownames(sample_info)))

phenoVec <- sample_info$ISRact
batchVec <- rownames(sample_info)


ads <- anota2seqDataSetFromMatrix(
  dataP = column_to_rownames(rowTumor_pol_gn, "geneID"),
  dataT = column_to_rownames(rowTumor_cyt_gn, "geneID"),
  phenoVec = phenoVec,
  batchVec = batchVec,
  dataType = "RNAseq",
  filterZeroGenes = TRUE, 
  normalize = TRUE,
  transformation = "TMM-log2",
  varCutOff = NULL)
