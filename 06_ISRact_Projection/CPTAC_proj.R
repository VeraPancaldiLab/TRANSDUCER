#Import library
library(tidyverse)
library(readxl)
library(biomaRt)
library(pdacmolgrad) #devtools::install_github("RemyNicolle/pdacmolgrad")
library(edgeR)
library(Hmisc)
#library(ggpubr)
################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

## CPTAC (Proteogenomics)
### upper quartile normalized and log2 data 
CPTAC_tumor_already_norm <- read_delim("data/PDAC_LinkedOmics_Data/mRNA_RSEM_UQ_log2_Tumor.cct", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE) %>%
  dplyr::rename(Gene = ...1)

## Raw counts obtained by email
CPTAC_tumor_raw <- read_delim("data/PDAC_LinkedOmics_Data/CPTAC-PDAC-raw-RSEM-expected-counts-gene-level.txt", 
                              delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(idx, matches("tumor")) %>%
  dplyr::rename_all(~str_remove(., '_tumor')) %>%
  dplyr::rename(Gene = idx)


### Metadata
low_purity_samples <- read_delim("data/PDAC_LinkedOmics_Data/low_purity_samples.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

clinical_data <- read_excel("data/PDAC_LinkedOmics_Data/mmc1.xlsx", 
                            sheet = "Clinical_data") %>%
  mutate(follow_up_days = as.numeric(follow_up_days)) %>%
  mutate(status = ifelse(vital_status == "Deceased", 2, 1))


Molecular_phenotype_data <- read_excel("data/PDAC_LinkedOmics_Data/mmc1.xlsx", 
                                       sheet = "Molecular_phenotype_data") %>% 
  mutate_at(vars(immune_deconv:`necrosis_(%OF_TUMOR_WITH_NECROSIS)_histology_estimate`, KRAS_VAF), as.numeric) %>%
  mutate_at(vars(Bailey:Moffitt), as.factor) %>% #change type to avoid errors
  dplyr::rename(sample = case_id)

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
norm_method = "upperquartile" #TMM | upperquartile | upperquartile_ogdata
#Variables for filter function
filter = FALSE # exclude just the low purity if FALSE
estimation_method = "histology" #histology | deconvolution
neoplastic_min = 0.2
acinar_max = 0.05
islet_max = 1
tumor_tissue_min = 0.8
################################################################################
# Choose raw / already normalized data
if (norm_method == "upperquartile_ogdata"){
  CPTAC_tumor__ <- CPTAC_tumor_already_norm
} else {
  CPTAC_tumor__ <- CPTAC_tumor_raw
}
  

# Translate EnsemblID to gene names
#Version 95 for Proteogenimics data
ensembl95 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 95)#listAttributes(ensembl95, page="feature_page")

annot_ensembl95 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_name'), mart = ensembl95)


#Add a Ensemble ID column to CPTAC_tumor__ (Ensembl release 95 according to the paper)
#https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=genes&hgta_track=mane&hgta_table=mane&hgta_doSchema=describe+table+schema

reverse_translate = deframe(annot_ensembl95[c("external_gene_name", "ensembl_gene_id")])

# Preprocessing
## filtering
if (filter == TRUE){
  CPTAC_tumor_ <- dataFilter(CPTAC_tumor__, clinical_data, Molecular_phenotype_data, estimation_method, neoplastic_min, acinar_max, islet_max, tumor_tissue_min)
} else {
  CPTAC_tumor_ <- dplyr::select(CPTAC_tumor__, -low_purity_samples$case_id)
}

## normalization and filtering
if (norm_method == "upperquartile_ogdata"){
  CPTAC_tumor <- CPTAC_tumor_
} else {
  CPTAC_tumor <- column_to_rownames(CPTAC_tumor_, "Gene") %>% DGEList() %>%
    calcNormFactors(method= norm_method) %>%
    cpm(log=TRUE) %>% as_tibble(rownames = "Gene")
}

## translate gene names
CPTAC_tumor$EnsemblID <- CPTAC_tumor$Gene %>%
  reverse_translate[.] %>%
  make.names(unique = TRUE)


## Create sample_info_CPTAC
type_pamg <- projectMolGrad(newexp = column_to_rownames(dplyr::select(CPTAC_tumor,-EnsemblID), "Gene"),  geneSymbols = CPTAC_tumor$Gene) %>%
  as_tibble(rownames = "sample")

sample_info_CPTAC <- dplyr::select(type_pamg, sample, ICGCrnaseq) %>% 
  dplyr::rename(PAMG = ICGCrnaseq) %>% 
  left_join(dplyr::select(Molecular_phenotype_data, sample, KRAS_VAF), "sample")

## Projecting datasets on this PCA
### CPTAC
#### transpose human data for projection
human_data <- CPTAC_tumor %>%
  pivot_longer(cols = 2:(length(CPTAC_tumor)-1), names_to = "case_id", values_to = "Expression") %>% 
  dplyr::filter(!str_detect(EnsemblID, 'NA')) %>% #remove eventual NA for gene Ensembl
  dplyr::select(-Gene) %>% 
  pivot_wider(names_from = EnsemblID, values_from = Expression, names_repair = "minimal") %>%
  column_to_rownames("case_id") 

### remove missing genes from the PCA before projection
filter_pca <- function(.data, objective){
  as_tibble(.data, rownames = "tmp") %>% 
    dplyr::filter(tmp %in% objective) %>% 
    column_to_rownames("tmp") %>%
    data.matrix()
}
pca_pdx$rotation <- filter_pca(pca_pdx$rotation, names(human_data)[-1])
pca_pdx$center <- filter_pca(pca_pdx$center, names(human_data)[-1])
pca_pdx$scale <- filter_pca(pca_pdx$scale, names(human_data)[-1])

projection_CPTAC <-  predict(pca_pdx, human_data) %>% 
  as_tibble(rownames = "sample")


CPTAC_PC1 <- arrange(projection_CPTAC, PC1) %>%
  mutate(PC1status = if_else(PC1 < quantile(projection_CPTAC$PC1, probs = 0.3333), "low_PC1",
                             if_else(PC1 < quantile(projection_CPTAC$PC1, probs = 0.6666), "medium_PC1", "high_PC1"))) %>%
  dplyr::select(sample, PC1,  PC1status) %>%
  left_join(sample_info_CPTAC, by = "sample")

#Plot pca and projections
#Add ISR status to PCA df
pca_full_df <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  inner_join(top_samples[,c("sample","ISRact")])

#-------------------------------------------------------------------------------
# Plot projecion PacaOmics and then projection Proteogenomics

show_projection <- bind_rows(as_tibble(pca_full_df) %>% mutate(dataset = "Sauyeun PDX"),
                             as_tibble(projection_CPTAC) %>% mutate(ISRact = "Unknown",
                                                                    dataset = "CPTAC"))

dplyr::mutate(show_projection, ISRact = str_replace(ISRact, 'ICA3', 'ISRact')) %>%
  dplyr::filter(dataset %in% c("Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
                ISRact %in% c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) %>% 
  ggplot(aes(x=PC1, y=PC2, color = ISRact, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Sauyeun PDX", "PACAOMICS PDX", "CPTAC", "CCLE")) +
  scale_color_discrete(limits = c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) + #c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue))
  ylim(-250,15)
#-------------------------------------------------------------------------------

# Plot comparisons with Basal/Classical and ISRact
## PAMG vs PC1
correlation_plotter(data = CPTAC_PC1, col1 = "PAMG", col2 = "PC1", data_name = "CPTAC")
