#Import library
library(tidyverse)
library(readxl)
library(biomaRt)
library(pdacmolgrad) #devtools::install_github("RemyNicolle/pdacmolgrad")
library(edgeR)
library(Hmisc)
library(ggrepel)
library(DEqMS)
library(matrixStats)
library(msigdbr)
library(fgsea)
library(ggpubr)
library(data.table)
library(ggplot2)
library(survival)
library(survminer)
library("AnnotationDbi")
library("illuminaHumanv4.db") 

################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

# Data loading
## Data
ICGC <- read_rds("data/ICGC/ICGC_AU.rds") %>% as_tibble(rownames = "Gene")
  
### Metadata
ICGC_clinical <- read_rds("data/ICGC/ICGC_AU_survival.rds")

### inherited sample info and top samples from Sauyeun_PDX
sample_info <- read_delim("data/Sauyeun_PDX/sample_info.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

top_samples <- arrange(sample_info, ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low_ICA3", "high_ICA3")) %>%
  arrange(sample)

# Import the projection
## ISRactPCA
pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")

################################################################################
# PARAMETERS
################################################################################

# Translate EnsemblID to gene names
probe_info <- AnnotationDbi::select(illuminaHumanv4.db,
                                    ICGC$Gene,
                                    c("SYMBOL", "ENSEMBL", "PROBEID"),
                                    keytype = "SYMBOL") %>%
  dplyr::rename(Gene = SYMBOL,
                EnsemblID = ENSEMBL)


# load annotation with Biomart
# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# 
# #listAttributes(ensembl75, page="feature_page")
# annot_ensembl <- getBM(attributes = c('ensembl_gene_id',
#                                         'external_gene_name'), mart = ensembl)
# 
# gene_to_ensembl = deframe(annot_ensembl[c( "external_gene_name", "ensembl_gene_id")])
# ensembl_to_gene = setNames(names(gene_to_ensembl), gene_to_ensembl)

## Deal with Duplicated EnsemblIDs and Gene names 
ICGC_ensembl <- left_join(ICGC, probe_info,by = "Gene") %>%
  dplyr::select(-c("Gene", "PROBEID")) %>%
  group_by(EnsemblID) %>%
  summarise_all(mean)

ICGC_gene <- ICGC %>%
  group_by(Gene) %>%
  summarise_all(mean)

# Calculate PAMG
type_pamg <- projectMolGrad(newexp = column_to_rownames(ICGC_gene, "Gene"),  geneSymbols = ICGC_gene$Gene) %>%
  as_tibble(rownames = "sample")

sample_info_Puleo <- dplyr::select(type_pamg, sample, ICGCarray) %>% 
  dplyr::rename(PAMG = ICGCarray) %>% 
  left_join(rownames_to_column(ICGC_clinical, "sample"), "sample")
