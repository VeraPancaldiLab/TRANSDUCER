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

################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

## CPTAC (Proteogenomics)
### upper quartile normalized and log2 data 
Puleo_already_norm2 <- read_rds("data/Puleo/Puleo2.rds") %>%
  as_tibble(rownames = "Gene") %>%
  mutate(Gene = strsplit(as.character(Gene), " /// ")) %>%
  unnest(Gene) %>%
  filter(Gene != "", 
         Gene != "---")


### Metadata
Survival_availability <- NULL
clinical_data <- read_rds("data/Puleo/Puleo_survival.rds")


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
norm_method = "upperquartile" #TMM | upperquartile | upperquartile_ogdata
