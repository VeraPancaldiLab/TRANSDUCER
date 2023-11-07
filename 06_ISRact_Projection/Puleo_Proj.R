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
library("hgu219.db")  

################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

# Data loading
## Data
Puleo_Shiny_Jeromme <- read_rds("data/Puleo/Puleo2.rds") %>%
  as_tibble(rownames = "Gene") %>%
  mutate(Gene = strsplit(as.character(Gene), " /// ")) %>%
  unnest(Gene) %>%
  filter(Gene != "", 
         Gene != "---")

## Translate EnsemblID to gene names
probe_info <- AnnotationDbi::select(hgu219.db,
                                    Puleo_Shiny_Jeromme$Gene,
                                    c("SYMBOL", "ENSEMBL"),
                                    keytype = "SYMBOL") %>%
  dplyr::rename(Gene = SYMBOL,
                EnsemblID = ENSEMBL)

## Deal with Duplicated EnsemblIDs and Gene names 
Puleo_Shiny_Jeromme_ensembl <- left_join(Puleo_Shiny_Jeromme, probe_info,by = "Gene") %>%
  dplyr::select(-"Gene") %>%
  group_by(EnsemblID) %>%
  summarise_all(mean)

Puleo_Shiny_Jeromme_gene <- Puleo_Shiny_Jeromme %>%
  group_by(Gene) %>%
  summarise_all(mean)

### Metadata
Survival_availability <- NULL
clinical_data <- read_rds("data/Puleo/Puleo_survival_repair.rds") %>%
  dplyr::rename(sample = ID)


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
# Choose normalization method (to be implemented?)
Puleo_gene <- Puleo_Shiny_Jeromme_gene
Puleo_ensembl <- Puleo_Shiny_Jeromme_ensembl

# Calculate PAMG
type_pamg <- projectMolGrad(newexp = column_to_rownames(Puleo_gene, "Gene"),  geneSymbols = Puleo_gene$Gene) %>%
  as_tibble(rownames = "sample")

sample_info_Puleo <- dplyr::select(type_pamg, sample, Puleo) %>% 
  dplyr::rename(PAMG = Puleo) %>% 
  left_join(clinical_data, "sample")