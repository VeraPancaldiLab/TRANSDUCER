#Import library
library(tidyverse)
library(biomaRt)
library(pdacmolgrad) #devtools::install_github("RemyNicolle/pdacmolgrad")
library(edgeR)
library(Hmisc)
#library(ggpubr)
################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

# Import datasets
## PACAOMICS
### Remy Nicolle's 2017 PDX data
RN2017_raw <- read_delim("data/PACAOMICs/Human-Tumor_rawcount_Transcriptome.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)
### Extended 73 sample cohort
load("data/PACAOMICs/Alexias_53PDX/PDX_HUMAN_RAW.RData")
PACAOMICs_90_raw <- as_tibble(x, rownames ="EnsemblID")

### inherited sample info and top samples from Sauyeun_PDX
sample_info <- read_delim("data/Sauyeun_PDX/sample_info.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

top_samples <-arrange(sample_info, ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low_ICA3", "high_ICA3")) %>%
  arrange(sample)


# Import the projection
## ISRactPCA
pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")
################################################################################
# PARAMETERS
PACAOMICS_raw = RN2017_raw #RN2017_raw | PACAOMICs_90_raw
norm_method = "upperquartile" #TMM | upperquartile
################################################################################

# Translate EnsemblID to gene names
## Version 75 for PDX data
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)#listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)

translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

# Processing and Normalization
## PACAOMICS
PACAOMICS_norm_ <- column_to_rownames(PACAOMICS_raw, "EnsemblID") %>%
  DGEList() %>%
  calcNormFactors(method= norm_method) %>%
  cpm(log=TRUE)

PACAOMICS_norm <- t(PACAOMICS_norm_) %>%
  as_tibble(rownames = "sample")

### Create sample_info_PACAOMICS
rownames(PACAOMICS_norm_) <- translate[rownames(PACAOMICS_norm_)]

type_pamg <- projectMolGrad(newexp = PACAOMICS_norm_,  geneSymbols =rownames(PACAOMICS_norm_)) %>%
  as_tibble(rownames = "sample")

sample_info_PACAOMICS <- dplyr::select(sample_info, -PAMG) %>%
  right_join(dplyr::select(type_pamg, sample, PDX)) %>% 
  dplyr::rename(PAMG = PDX)

# Projection
filter_pca <- function(.data, objective){
  as_tibble(.data, rownames = "tmp") %>% 
    dplyr::filter(tmp %in% objective) %>% 
    column_to_rownames("tmp") %>%
    data.matrix()
}

pca_pdx$rotation <- filter_pca(pca_pdx$rotation, names(PACAOMICS_norm)[-1])
pca_pdx$center <- filter_pca(pca_pdx$center, names(PACAOMICS_norm)[-1])
pca_pdx$scale <- filter_pca(pca_pdx$scale, names(PACAOMICS_norm)[-1])


projection_PACAOMICS <- predict(pca_pdx, PACAOMICS_norm) %>%
  as_tibble() %>%
  mutate(sample = PACAOMICS_norm$sample, .before = 1) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace(ISRact, is.na(ISRact) & sample %in% sample_info$sample, "medium_ICA3"),
         ISRact = replace_na(ISRact, "Unknown"))

PACAOMICS_PC1 <- arrange(projection_PACAOMICS, PC1) %>% #! wrong! remove the last part of medium or the full object
  dplyr::filter(!sample=="!") %>%
  mutate(PC1status = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3/3))), labels = c("low_PC1", "medium_PC1", "high_PC1"), include.lowest = TRUE)) %>%
  dplyr::select(sample, PC1,  PC1status) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace_na(ISRact, "medium_ICA3")) %>%
  inner_join(sample_info_PACAOMICS[, c("sample","PAMG", "ICA3")], by = "sample")



# Plot pca and projections
## Add ISR status to PCA df
pca_full_df <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  inner_join(top_samples[,c("sample","ISRact")])

#-------------------------------------------------------------------------------
# Plot projecion

show_projection <- bind_rows(as_tibble(pca_full_df) %>% mutate(dataset = "Sauyeun PDX"),
                             mutate(projection_PACAOMICS, sample = paste0(sample, "_OG"),
                                    dataset = "PACAOMICS PDX"))

dplyr::mutate(show_projection, ISRact = str_replace(ISRact, 'ICA3', 'ISRact')) %>%
  dplyr::filter(dataset %in% c("Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
                ISRact %in% c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) %>% 
  ggplot(aes(x=PC1, y=PC2, color = ISRact, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Sauyeun PDX", "PACAOMICS PDX", "CPTAC", "CCLE")) +
  scale_color_discrete(limits = c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) #c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue)
#-------------------------------------------------------------------------------

# Plot comparisons with Basal/Classical and ISRact
## ISR vs PC1
correlation_plotter(data = PACAOMICS_PC1, col1 = "ICA3", col2 = "PC1", data_name = "PACAOMICs PDX")
## PAMG vs PC1
correlation_plotter(data = PACAOMICS_PC1, col1 = "PAMG", col2 = "PC1", data_name = "PACAOMICs PDX")
## ISR vs PAMG
correlation_plotter(data = PACAOMICS_PC1, col1 = "ICA3", col2 = "PAMG", data_name = "PACAOMICs PDX")
