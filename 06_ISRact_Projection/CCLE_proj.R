#Import library
library(tidyverse)
library(biomaRt)
library(pdacmolgrad) #devtools::install_github("RemyNicolle/pdacmolgrad")
library(edgeR)
library(Hmisc)
library(ggrepel)
################################################################################
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

# Import datasets
## CCLE 
ccle_raw <- read_delim(file="data/CCLE/CCLE_RNAseq_genes_counts_20180929.gct", skip=2) %>%
  dplyr::rename(EnsemblID = Name, GeneName = Description) %>% dplyr::mutate(EnsemblID = str_remove(EnsemblID, '\\.[0-9]*$'))

ccle_info <- read_csv("data/CCLE/primary-screen-cell-line-info.csv") %>%
  dplyr::mutate(ISRact = if_else(ccle_name %in% c("ASPC1_PANCREAS", "PATU8902_PANCREAS", "PATU8988T_PANCREAS", "MIAPACA2_PANCREAS"),
                                 if_else(ccle_name %in% c("ASPC1_PANCREAS", "PATU8902_PANCREAS"),
                                         "high_ISRact", "low_ISRact"),
                                 NA))

### Sample info and top samples from Sauyeun_PDX
Sauyeun_info <- read_delim("data/Sauyeun_PDX/sample_info.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

top_samples <- arrange(Sauyeun_info, ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low_ICA3", "high_ICA3")) %>%
  arrange(sample)

# Import the projection
## ISRactPCA
pca_pdx <- read_rds("data/Classifiers/pca_pdx_ENZO.RDS")

################################################################################
# PARAMETERS
keep_CCLE = "PANCREAS" #ALL | PANCREAS
norm_method = "upperquartile" #TMM | upperquartile
################################################################################

# Processing and Normalization
## CCLE
if (keep_CCLE == "PANCREAS"){
  keep <- dplyr::filter(ccle_info, primary_tissue == "pancreas") %>% 
    dplyr::select(ccle_name) %>%
    deframe()
  ccle_raw <- dplyr::select(ccle_raw, EnsemblID, GeneName, all_of(keep))
}
ccle_norm_ <- dplyr::select(ccle_raw, -GeneName) %>%
  column_to_rownames("EnsemblID") %>%
  DGEList() %>%
  calcNormFactors(method= norm_method) %>%
  cpm(log=TRUE) 

ccle_norm <- ccle_norm_ %>%
  t() %>%
  as_tibble(rownames = "ccle_name")

### create version with Gene Names
ccle_translate = deframe(ccle_raw[c("EnsemblID", "GeneName")])

ccle_norm_gn <- as_tibble(ccle_norm_, rownames = "EnsemblID") %>%
  mutate(GeneName = ccle_translate[EnsemblID]) %>% 
  dplyr::select(-EnsemblID) %>% group_by(GeneName) %>% 
  summarise_all(sum)

### Add PAMG status
type_pamg <- projectMolGrad(newexp = column_to_rownames(ccle_norm_gn, "GeneName"),  geneSymbols = ccle_norm_gn$GeneName) %>%
  as_tibble(rownames = "ccle_name")

ccle_info <- dplyr::select(type_pamg, ccle_name, PDX) %>% 
  right_join(ccle_info, "ccle_name") %>%
  dplyr::rename(PAMG = PDX)


# Projection
## CCLE
ccle_norm_minval <- min(ccle_norm_) # same minimum value (unlike TMM?)
ccle_missing_genes <- rownames(pca_pdx$rotation)[(rownames(pca_pdx$rotation) %in% names(ccle_norm)) == F]

ccle_missing_data <- as_tibble(matrix(ccle_norm_minval,
                                      ncol = length(ccle_missing_genes),
                                      nrow = nrow(ccle_norm),
                                      dimnames = list(ccle_norm$ccle_name, ccle_missing_genes)),
                               rownames = "ccle_name")

projection_ccle <- predict(pca_pdx, inner_join(ccle_norm, ccle_missing_data, by="ccle_name")) %>% 
  as_tibble() %>%
  mutate(ccle_name = ccle_norm$ccle_name, .before = 1) %>%
  left_join(ccle_info[,c("ccle_name","primary_tissue", "ISRact", "PAMG")], by="ccle_name")

ccle_PC1 <- arrange(projection_ccle, PC1) %>% 
  mutate(PC1status = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3/3))), labels = c("low_PC1", "medium_PC1", "high_PC1"), include.lowest = TRUE)) %>%
  dplyr::select(ccle_name, PC1,  PC1status)


# Plot pca and projections
## Add ISR status to PCA df
pca_full_df <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  inner_join(top_samples[,c("sample","ISRact")])

#-------------------------------------------------------------------------------
# Plot projecion

show_projection <- bind_rows(as_tibble(pca_full_df) %>% mutate(dataset = "Sauyeun PDX"),
                             as_tibble(projection_ccle %>% mutate(ISRact = if_else(is.na(ISRact), "Unknown", ISRact),
                                                                  dataset = "CCLE")))

dplyr::mutate(show_projection, ISRact = str_replace(ISRact, 'ICA3', 'ISRact')) %>%
  dplyr::filter(dataset %in% c("Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
                ISRact %in% c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) %>% 
  ggplot(aes(x=PC1, y=PC2, color = ISRact, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Sauyeun PDX", "PACAOMICS PDX", "CPTAC", "CCLE")) +
  scale_color_discrete(limits = c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) #c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue)
#-------------------------------------------------------------------------------

# Plot comparisons with Basal/Classical and ISRact
## PAMG vs PC1
correlation_plotter(data = projection_ccle, col1 = "PAMG", col2 = "PC1", data_name = "CCLE")


#-------------------------------------------------------------------------------
# CCLE plot with names and colored according to any gene expression
pivot_longer(ccle_norm_gn,-GeneName, names_to = "ccle_name", values_to = "expression") %>%
  pivot_wider(id_cols = "ccle_name", names_from = "GeneName", values_from = "expression") %>%
  inner_join(projection_ccle, by="ccle_name") %>% 
  ggplot(aes(x=PC1, y=GATA6, label = str_remove(ccle_name, "_PANCREAS"))) +
  geom_point(aes(color = ISRact)) + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel()
#-------------------------------------------------------------------------------
