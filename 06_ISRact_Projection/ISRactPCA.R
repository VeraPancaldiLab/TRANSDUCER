#Import library
library(tidyverse)
library(biomaRt)
library(edgeR)
library(factoextra)
library(ggpubr)
################################################################################
#List the n most absolute correlated genes with ICA3
#'@description
#'this function calculate the chosen absolute correlation between genes and ICA3 component 
#'with either full samples (best overall) or ISR samples (best subset) and then subset the top n
#'@param normdf data frame of normalized RNA-seq count data. Lines are samples and columns are genes
#'@param select_method "bestsubset" to select ISR samples before doing the correlation or 
#'"bestoverall" to keep all the samples
#'@param n number of gene to keep
#'@param cormethod "pearson" or "spearman"
#'@param top_samples data frame of top sample info
#'@param sample_info data frame of every sample info

mostCorrGenes <- function(normdf, select_method, n, cormethod, top_samples, sample_info){
  
  #Keep only top samples if selection method is bestsubset
  if (select_method == "bestsubset"){
    normdf <- dplyr::filter(normdf, sample %in% top_samples$sample) %>%
      arrange(sample) 
    sample_df <- top_samples
  } else if (select_method == "bestoverall"){
    normdf <- arrange(normdf, sample) 
    sample_df <- sample_info
  }
  
  #Calculate absolute correlation for each gene
  genes_ic3_cor <- normdf %>% 
    dplyr::filter(!sample=="!") %>%
    column_to_rownames("sample") %>% #select continuous variables
    as.matrix() %>%
    cor(y = sample_df$ICA3, method = cormethod) %>%
    as.data.frame() 
  
  #keep the 1000 most absolute correlated genes
  top_genes <- mutate(genes_ic3_cor, V1 = abs(V1)) %>%
    arrange(desc(V1)) %>%
    dplyr::slice(1:n) %>%
    rownames_to_column(var = "EnsemblID") 
  
  
  return(top_genes)
}
################################################################################


##devtools::install_github("RemyNicolle/pdacmolgrad")
setwd("~/Documents/02_TRANSDUCER/06_ISRact_Projection/")

#Import filter function
source("src/human_cohort_data_filter.R")
source("src/correlation_plotter.R")

#Import datasets
## Sauyeun PDX
### metadata
sample_info <- read_delim("data/Sauyeun_PDX/sample_info.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

sample_names_equivalences <- read_delim("data/Sauyeun_PDX/sample_names_equivalences.tsv", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)

### raw count processed by Jacobo
rawTumor_Cyt <- read_delim("data/Sauyeun_PDX/rawTumor_Cyt.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE) %>%
  dplyr::select(-`!`)

### raw count processed by Remy Nicolle
geneCount_raw_28s_totalRNA <- read_delim("data/Sauyeun_PDX/geneCount_raw_28s_totalRNA.tsv", 
                                         delim = "\t", escape_double = FALSE, 
                                         trim_ws = TRUE)  %>%
  rename_with(~coalesce(sample_names_equivalences$CITID[match(., sample_names_equivalences$fastq_name)], .)) %>%
  dplyr::select(-`!`)

### best subset selected data direct input on PCA
ISRact_data_bestsubset <- read_delim("data/Sauyeun_PDX/5vs5_training_TMM_nonstandardized_bestsubset.tsv", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)  %>%
  dplyr::rename(sample = rowname)

### best overall selected data direct input on PCA
ISRact_data_bestoverall <- read_delim("data/Sauyeun_PDX/5vs5_training_TMM_nonstandardized_bestoverall.tsv", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE) %>%
  dplyr::rename(sample = rowname)


################################################################################
# PARAMETERS
Sauyeun_raw = rawTumor_Cyt #geneCount_raw_28s_totalRNA | rawTumor_Cyt
cormethod = "pearson"
keepgenes = 1000
select_method = "bestsubset" #bestsubset | bestoverall
norm_method = "upperquartile" #TMM | upperquartile

#Variables for filter function
filter = FALSE # exclude just the low purity if FALSE
estimation_method = "histology" #histology | deconvolution
neoplastic_min = 0.2
acinar_max = 0.05
islet_max = 1
tumor_tissue_min = 0.8
################################################################################

# Translate EnsemblID to gene names
## Version 75 for PDX data
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)#listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)

translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

#Add a Gene names column to Sauyeun_raw
Sauyeun_raw$EnsemblID %>% 
  translate[.] %>%
  make.names(unique = TRUE) -> Sauyeun_raw$Genenames

# Processing and Normalization
## Normalize raw counts
Sauyeun_norm_ <- dplyr::select(Sauyeun_raw, -Genenames) %>%
  column_to_rownames( "EnsemblID") %>%
  DGEList() %>%
  calcNormFactors(method= norm_method) %>%
  cpm(log=TRUE) 

Sauyeun_norm <- Sauyeun_norm_ %>%
  t() %>%
  as_tibble(rownames = "sample")

### Visualize it
pivot_longer(Sauyeun_norm,
             cols = 2:length(Sauyeun_norm),
             names_to = "EnsemblID",
             values_to = "Expression", 
             names_repair = "unique") %>%
  ggplot(data , mapping = aes(x = sample, y = Expression, fill = sample))+
  geom_violin() +
  theme(legend.position = "none") +
  coord_flip()

## sample and gene filtering
### get 5 and 5 most extreme samples for ICA3 
top_samples <-arrange(sample_info, ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low_ICA3", "high_ICA3")) %>%
  arrange(sample)

ggdensity(
  sample_info, x = "ICA3",
  rug = TRUE
) + geom_rug(data = top_samples, 
             inherit.aes = F, 
             aes(x = ICA3, color = ISRact))

### get the most correlated genes with ICA3
top_genes <- mostCorrGenes(Sauyeun_norm,
                           select_method,
                           keepgenes,
                           cormethod,
                           top_samples,
                           sample_info)

### Final subset
Sauyeun_norm_subset <- dplyr::select(Sauyeun_norm, sample, top_genes$EnsemblID) %>% 
  dplyr::filter(sample %in% top_samples$sample) %>%
  arrange(sample)


# PCA 
pca_pdx <- Sauyeun_norm_subset %>% 
  column_to_rownames( "sample") %>% 
  scale() %>%
  prcomp()


## Visualize 
fviz_eig(pca_pdx, addlabels = TRUE, ylim = c(0, 50)) #visualize explained variance

fviz_pca_ind(pca_pdx,
             col.ind = top_samples$ISRact, # color by groups
             palette = c("#00FF00", "#FF0000"),
             #addEllipses = TRUE, # Ellipses 
             legend.title = "Groups",
             invisible="quali"
)

## export PCA for projection
write_rds(pca_pdx, "data/Classifiers/pca_pdx.RDS")

## create object for further comparisons
projection_Sauyeun <- predict(pca_pdx, Sauyeun_norm) %>%
  as_tibble() %>%
  mutate(sample = Sauyeun_norm$sample, .before = 1) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace_na(ISRact, "medium_ICA3")) 

Sauyeun_PC1 <- arrange(projection_Sauyeun, PC1) %>% 
  dplyr::filter(!sample=="!") %>%
  mutate(PC1status = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3/3))), labels = c("low_PC1", "medium_PC1", "high_PC1"), include.lowest = TRUE)) %>%
  dplyr::select(sample, PC1,  PC1status) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace_na(ISRact, "medium_ICA3")) %>%
  inner_join(sample_info[, c("sample","PAMG", "ICA3", "Diabetes")], by = "sample")

# Plot comparisons with Basal/Classical and ISRact
## ISR vs PC1
correlation_plotter(data = Sauyeun_PC1, col1 = "ICA3", col2 = "PC1", data_name = "Sauyeun PDX")
## PAMG vs PC1
correlation_plotter(data = Sauyeun_PC1, col1 = "PAMG", col2 = "PC1", data_name = "Sauyeun PDX")
## ISR vs PAMG
correlation_plotter(data = Sauyeun_PC1, col1 = "ICA3", col2 = "PAMG", data_name = "Sauyeun PDX")
## ISR vs IFNsign
correlation_plotter(data = Sauyeun_PC1, col1 = "ICA3", col2 = "IFNsign", data_name = "Sauyeun PDX")


