# ISRactPCA projection into Maurer 2019 LMD cohort
################################################################################
#' Import different files
#' Filter data and project ISRactPCA
#' Produce projection scatterplot
#' Compare ISRactPCA with PAMG score
#' For more information run 
################################################################################

# Import libraries
library(tidyverse)
library(biomaRt)
library(pdacmolgrad) # devtools::install_github("RemyNicolle/pdacmolgrad")
library(edgeR)
library(Hmisc)
library(ggrepel)
################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

# Import datasets
## CPTAC (Proteogenomics)
### Load data
load("data/Maurer2019/Start Maurer EvsSvsB.RData")

## ISRactPCA projection
pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")
### inherited sample info and top samples from Shin PDX training data
sample_info <- read_delim("data/Sauyeun_PDX/sample_info.tsv",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE
)

top_samples <- arrange(sample_info, ICA3) %>%
  dplyr::slice(unique(c(1:5, n() - 0:4))) %>%
  mutate(ISRact_bin = ifelse(ICA3 < 0, "low_ISRact", "high_ISRact")) %>%
  dplyr::rename(ISRact = "ICA3") %>%
  arrange(sample)

################################################################################
# PARAMETERS
norm_method <- "upperquartile" # TMM | upperquartile
subset_celltype <- "epithelium" # epithelium | stroma | bulk
################################################################################

# Translate EnsemblID to gene names
## Version 109
ensembl75 <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  version = 109,
  host = "https://grch37.ensembl.org"
) # listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c(
  "ensembl_gene_id",
  "hgnc_symbol"
), mart = ensembl75)

translate <- deframe(annot_ensembl75[c("ensembl_gene_id", "hgnc_symbol")])
inverse_translate <- deframe(annot_ensembl75[c("hgnc_symbol", "ensembl_gene_id")])

# Preprocessing
## Subset cell type in Maurer dataset
sample_type <- dplyr::filter(phenoEvsS, Compartment == subset_celltype)
maurer_rawcount <- dplyr::select(raw_E_Vs_S, sample_type$Sample) %>%
  dplyr::mutate(EnsemblID = inverse_translate[rownames(.)]) %>%
  dplyr::group_by(EnsemblID) %>%
  dplyr::summarise_all(sum) %>%
  dplyr::filter(!is.na(EnsemblID)) %>%
  column_to_rownames("EnsemblID")


maurer_normcount_ <- DGEList(maurer_rawcount) %>%
  calcNormFactors(method = norm_method) %>%
  cpm(log = TRUE)

maurer_normcount <- t(maurer_normcount_) %>%
  as_tibble(rownames = "sample")


### Create sample_info_maurer
type_pamg <- projectMolGrad(newexp = maurer_normcount_, geneSymbols = translate[rownames(maurer_normcount_)]) %>%
  as_tibble(rownames = "sample")

maurer_info <- dplyr::select(type_pamg, sample, ICGCrnaseq) %>%
  dplyr::rename(PAMG = ICGCrnaseq)

# Projection
## Maurer
### Add missing genes to expression before the PCA
filter_pca <- function(.data, objective) {
  as_tibble(.data, rownames = "tmp") %>%
    dplyr::filter(tmp %in% objective) %>%
    column_to_rownames("tmp") %>%
    data.matrix()
}

pca_pdx$rotation <- filter_pca(pca_pdx$rotation, names(maurer_normcount)[-1])
pca_pdx$center <- filter_pca(pca_pdx$center, names(maurer_normcount)[-1])
pca_pdx$scale <- filter_pca(pca_pdx$scale, names(maurer_normcount)[-1])

Maurer_PCAspace <- predict(pca_pdx, column_to_rownames(maurer_normcount, "sample")) %>%
  as_tibble(rownames = "sample") %>%
  left_join(maurer_info[, c("sample", "PAMG")], by = "sample") %>%
  dplyr::rename(ISRactPCA = "PC1") %>%
  mutate(
    ISRactPCA_bin = if_else(ISRactPCA < quantile(ISRactPCA, probs = 0.3333), "low_ISRactPCA",
      if_else(ISRactPCA < quantile(ISRactPCA, probs = 0.6666), "intermediate_ISRactPCA", "high_ISRactPCA")
    ),
    ISRactPCA_bin = fct(ISRactPCA_bin, levels = c("low_ISRactPCA", "intermediate_ISRactPCA", "high_ISRactPCA"))
  )

Maurer_ISRact_projection <- arrange(Maurer_PCAspace, ISRactPCA) %>%
  dplyr::select(sample, ISRactPCA, ISRactPCA_bin)


# Plot pca and projections
## Add ISR status to PCA df
pca_training <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  dplyr::rename(ISRactPCA = "PC1") %>%
  inner_join(top_samples[, c("sample", "ISRact_bin")])

#-------------------------------------------------------------------------------
# Plot projecion

show_projection <- bind_rows(
  as_tibble(pca_training) %>% mutate(dataset = "Sauyeun PDX"),
  as_tibble(Maurer_PCAspace %>% mutate(
    ISRact_bin = "Unknown",
    dataset = "Maurer"
  ))
)

dplyr::filter(
  show_projection, dataset %in% c("Sauyeun PDX", "PACAOMICS PDX", "CPTAC", "CCLE", "Maurer"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
  ISRact_bin %in% c("low_ISRact", "high_ISRact", "intermediate_ISRact", "Unknown")
) %>%
  ggplot(aes(x = ISRactPCA, y = PC2, color = ISRact_bin, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Sauyeun PDX", "PACAOMICS PDX", "CPTAC", "CCLE", "Maurer")) +
  scale_color_discrete(limits = c("low_ISRact", "high_ISRact", "intermediate_ISRact", "Unknown")) # c('low_ISRact', 'high_ISRact', 'intermediate_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue)
#-------------------------------------------------------------------------------

# Plot comparisons with Basal/Classical and ISRact
## PAMG vs ISRactPCA
correlation_plotter(data = Maurer_PCAspace, col1 = "PAMG", col2 = "ISRactPCA", data_name = "Maurer")
