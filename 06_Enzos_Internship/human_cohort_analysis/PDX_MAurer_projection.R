#Maurer projection

################################################################################
##Upper-quantile and log2 normalize Sauyeun pdx data
#Redo the Jacobo's PCA with 10 most extreme 
#sample regarding ISR activation for PDX data
#Project Maurer cohort dataset in this PCA
################################################################################

#Import library
library(readr)
library(readxl)
library(biomaRt)
library(tidyverse)
library(Hmisc)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(edgeR)
library(rcompanion)
library(corrr)
library(rstatix)
library(survminer)
library(pdacmolgrad)


setwd("~/Documents/02_TRANSDUCER/06_Enzos_Internship/human_cohort_analysis/")

#Import dataset
#Jacobo's PDX RNA-seq raw count processed by himself
rawTumor_Cyt <- read_delim("Data/rawTumor_Cyt.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

#Maurer raw count and sample info
load("Data/Start Maurer EvsSvsB.RData")

#Information of Jacobo's PDX samples including ICA results
sample_info <- read_delim("Data/sample_info.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

################################################################################
#Variables for normalization and gene selection
rawdata = rawTumor_Cyt 
cormethod = "pearson"
keepgenes = 1000
select_method = "bestsubset" #bestsubset | bestoverall
norm_method = "upperquartile" #TMM | upperquartile
#Variable to filter Maurer data
subset_celltype = "epithelium" # epithelium | stroma | bulk
################################################################################


####Translate EnsemblID to gene names####
# load annotation with Biomart
#Version 75 for PDX data
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)#listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = ensembl75)


#Add a Gene names column to rawdata
translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

rawdata$EnsemblID %>% 
  translate[.] %>%
  make.names(unique = TRUE) -> rawdata$Genenames

#Subset cell type in Maurer dataset
sample_type <- dplyr::filter(phenoEvsS, Compartment == subset_celltype)

maurer_rawcount <- dplyr::select(raw_E_Vs_S, sample_type$Sample)


####Normalize raw counts####
normdata <- dplyr::select(rawdata, -EnsemblID) %>%
  column_to_rownames( "Genenames") %>%
  DGEList() %>%
  calcNormFactors(method= norm_method) %>%
  cpm(log=TRUE) %>%
  t() %>%
  as_tibble(rownames = "sample")

maurer_normcount_ <- DGEList(maurer_rawcount) %>%
  calcNormFactors(method= norm_method) %>%
  cpm(log=TRUE) 

maurer_normcount <- %>%
  t(maurer_normcount_) %>%
  as_tibble(rownames = "sample")

### Create sample_info_maurer

type_pamg <- projectMolGrad(newexp = maurer_normcount_,  geneSymbols = rownames(maurer_normcount_)) %>%
  as_tibble(rownames = "sample")
sample_info_Maurer <- dplyr::select(type_pamg, sample, ICGCrnaseq) %>% 
  dplyr::rename(PAMG = ICGCrnaseq)


####Subset top genes and extreme samples####
#List 10 most extreme samples for ICA3 
ICA3_density <- ggdensity(
  sample_info, x = "ICA3", 
  rug = TRUE
) 
ICA3_density

top_samples <- sample_info %>%
  dplyr::arrange(ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low_ICA3", "high_ICA3")) %>%
  dplyr::arrange(sample) 

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
    dplyr::arrange(desc(V1)) %>%
    dplyr::slice(1:n) %>%
    rownames_to_column(var = "Gene") 
  
  
  return(top_genes)
}

#Apply function to select top genes
top_genes <- mostCorrGenes(normdata, select_method, keepgenes, cormethod, top_samples, sample_info)

#Check if they are present in mRNA_tumor
table(top_genes$Gene %in% rownames(maurer_rawcount)) #only 692 common

#Final subset
normdata_topg_s <- dplyr::select(normdata, sample, top_genes$Gene) %>% #to keep Jacobo's genes: any_of(colnames(ISRact_data_bestsubset))
  dplyr::filter(sample %in% top_samples$sample) %>%
  arrange(sample)

#Run the PCA 
pca_pdx <- normdata_topg_s %>% 
  column_to_rownames( "sample") %>% 
  dplyr::select(any_of(rownames(maurer_rawcount))) %>%
  scale() %>%
  prcomp()


#Visualize dimensions explained variance
fviz_eig(pca_pdx, addlabels = TRUE, ylim = c(0, 50)) #visualize explained variance

#Visualize individuals
fviz_pca_ind(pca_pdx,
             col.ind = top_samples$ISRact, # color by groups
             palette = c("#00FF00", "#FF0000"),
             #addEllipses = TRUE, # Ellipses 
             legend.title = "Groups",
             invisible="quali"
)

##Projecting Maurer dataset on this PCA
projection_Maurer <- predict(pca_pdx, maurer_normcount) %>%
  as_tibble() %>%
  mutate(sample = maurer_normcount$sample, .before = 1) 

#Plot pca and projections
#Add ISR status to PCA df
pca_full_df <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  inner_join(top_samples[,c("sample","ISRact")])

#Plot
ggplot(projection_Maurer, aes(x=PC1, y=PC2)) +
  geom_point() +
  geom_point(data = pca_full_df, aes(color =  ISRact))

#write PC1 weight
Maurer_PC1 <- arrange(projection_Maurer, PC1) %>% 
  mutate(PC1status = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3/3))), labels = c("low_PC1", "medium_PC1", "high_PC1"), include.lowest = TRUE)) %>%
  dplyr::select(sample, PC1,  PC1status)  %>%
  left_join(sample_info_Maurer, by="sample")
  
write.csv(Maurer_PC1, "Maurer_bulk_PC1.csv", row.names = FALSE)

# check correlation with PAMG
## Function to enhance correlation visualization
correlation_plotter <- function(data = Sauyeun_PC1, col1 = "PC1", col2 = "PAMG", data_name = "sauyeun PDX"){
  corr_spearman <- rcorr(data[[col1]], data[[col2]], type = "spearman")
  corr_pearson <- rcorr(data[[col1]], data[[col2]], type = "pearson")
  stats <- paste0("Spearman: R = ", round(corr_spearman$r["x","y"], 2), ", pval = ", round(corr_spearman$P["x","y"], 4),
                  "\nPearson: R = ", round(corr_pearson$r["x","y"], 2), ", pval = ", round(corr_pearson$P["x","y"], 4))
  
  ggplot(data) +
    aes_string(col1, col2) +
    geom_point(shape = 16, size = 2, show.legend = FALSE) +
    geom_smooth(method=lm) +
    theme_minimal() +
    labs(title = paste0("Comparison between ", col1, " and ", col2, " in ", data_name),
         subtitle = stats)
}

## PAMG vs PC1
correlation_plotter(data = Maurer_PC1, col1 = "PAMG", col2 = "PC1", data_name = "Maurer epithelial")
