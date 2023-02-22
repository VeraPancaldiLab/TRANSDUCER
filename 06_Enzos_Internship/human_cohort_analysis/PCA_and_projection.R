### Internship mission 2 
#(old name is PDX_human_projection.R)
################################################################################

#Upper-quantile and log2 normalize pdx data
#Redo the Jacobo's PCA with 10 most extreme 
#sample regarding ISR activation for PDX data
#Try to project human cohort dataset in this PCA
#GSEA reactome PDX data

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
library(GSVA)
##devtools::install_github("RemyNicolle/pdacmolgrad")


setwd("~/Documents/02_TRANSDUCER/06_Enzos_Internship/human_cohort_analysis/Data/")

#Import filter function
source("../human_cohort_data_filter.R")

#Import datasets
## Sauyeun PDX
### metadata
sample_info <- read_delim("sample_info.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

sample_names_equivalences <- read_delim("sample_names_equivalences.tsv", 
                                        delim = "\t", escape_double = FALSE, 
                                        trim_ws = TRUE)

### raw count processed by Jacobo
rawTumor_Cyt <- read_delim("rawTumor_Cyt.tsv", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE) %>%
  dplyr::select(-`!`)

### raw count processed by Remy Nicolle
geneCount_raw_28s_totalRNA <- read_delim("geneCount_raw_28s_totalRNA.tsv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)  %>%
  rename_with(~coalesce(sample_names_equivalences$CITID[match(., sample_names_equivalences$fastq_name)], .)) %>%
  dplyr::select(-`!`)

### best subset selected data direct input on PCA
ISRact_data_bestsubset <- read_delim("5vs5_training_TMM_nonstandardized_bestsubset.tsv", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)  %>%
  dplyr::rename(sample = rowname)

### best overall selected data direct input on PCA
ISRact_data_bestoverall <- read_delim("5vs5_training_TMM_nonstandardized_bestoverall.tsv", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE) %>%
  dplyr::rename(sample = rowname)


## CPTAC (Proteogenomics)
### upper quartile normalized and log2 data 
CPTAC_tumor_raw <- read_delim("PDAC_LinkedOmics_Data/mRNA_RSEM_UQ_log2_Tumor.cct", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)%>%
     dplyr::rename(Gene = ...1)

### Metadata
low_purity_samples <- read_delim("low_purity_samples.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

clinical_data <- read_excel("mmc1.xlsx", 
                            sheet = "Clinical_data") %>%
  mutate(follow_up_days = as.numeric(follow_up_days)) %>%
  mutate(status = ifelse(vital_status == "Deceased", 2, 1))


Molecular_phenotype_data <- read_excel("mmc1.xlsx", 
                                       sheet = "Molecular_phenotype_data") %>% 
  mutate_at(vars(immune_deconv:`necrosis_(%OF_TUMOR_WITH_NECROSIS)_histology_estimate`, KRAS_VAF), as.numeric) %>%
  mutate_at(vars(Bailey:Moffitt), as.factor) #change type to avoid errors

## PACAOMICS
### Remy Nicolle's 2017 PDX data
RN2017_raw <- read_delim("Human-Tumor_rawcount_Transcriptome.tsv", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)
### Extended 73 sample cohort
load("Alexias_53PDX/PDX_HUMAN_RAW.RData")
PACAOMICs_90_raw <- as_tibble(x, rownames ="EnsemblID")

## CCLE 
ccle_raw <- read_delim(file="CCLE/CCLE_RNAseq_genes_counts_20180929.gct", skip=2) %>%
  dplyr::rename(EnsemblID = Name, GeneName = Description) %>% dplyr::mutate(EnsemblID = str_remove(EnsemblID, '\\.[0-9]*$'))

ccle_info <- read_csv("CCLE/primary-screen-cell-line-info.csv") %>%
  dplyr::mutate(ISRact = if_else(ccle_name %in% c("ASPC1_PANCREAS", "PATU8902_PANCREAS", "PATU8988T_PANCREAS", "MIAPACA2_PANCREAS"),
                                 if_else(ccle_name %in% c("ASPC1_PANCREAS", "PATU8902_PANCREAS"),
                                         "high_ISRact", "low_ISRact"),
                                 NA))

## Other data
### Espinet et al 2019 IFNsign
IFNsign_geneset <- read_tsv("IFNsign_Espinet.tsv",col_names = "IFNsign")

################################################################################
# PARAMETERS
Sauyeun_raw = rawTumor_Cyt #geneCount_raw_28s_totalRNA | rawTumor_Cyt
PACAOMICS_raw = PACAOMICs_90_raw #RN2017_raw | PACAOMICs_90_raw
keep_CCLE = "PANCREAS" #ALL | PANCREAS 
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

#Version 95 for Proteogenimics data
ensembl95 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 95)#listAttributes(ensembl95, page="feature_page")

annot_ensembl95 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_name'), mart = ensembl95)

#Add a Gene names column to Sauyeun_raw
translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

Sauyeun_raw$EnsemblID %>% 
  translate[.] %>%
  make.names(unique = TRUE) -> Sauyeun_raw$Genenames


#Add a Ensemble ID column to CPTAC_tumor_raw (Ensembl release 95 according to the paper)
#https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=genes&hgta_track=mane&hgta_table=mane&hgta_doSchema=describe+table+schema

reverse_translate = deframe(annot_ensembl95[c("external_gene_name", "ensembl_gene_id")])

CPTAC_tumor_raw$EnsemblID <- CPTAC_tumor_raw$Gene %>% 
  reverse_translate[.] %>%
  make.names(unique = TRUE)

# Processing and Normalization
## Sauyeun PDX
### Normalize raw counts
Sauyeun_norm_ <- dplyr::select(Sauyeun_raw, -Genenames) %>%
  column_to_rownames( "EnsemblID") %>%
  DGEList() %>%
  calcNormFactors(method= norm_method) %>%
  cpm(log=TRUE) 

Sauyeun_norm <- Sauyeun_norm_ %>%
  t() %>%
  as_tibble(rownames = "sample")

### Visualize it
data <- pivot_longer(Sauyeun_norm, cols = 2:length(Sauyeun_norm), names_to = "EnsemblID", values_to = "Expression", names_repair = "unique")

g1 <- ggplot(data , mapping = aes(x = sample, y = Expression, fill = sample))+
  geom_violin()+
  theme(legend.position = "none")+
  coord_flip()

g1

### Keep common genes with RN2017_raw 
Sauyeun_norm <- dplyr::select(Sauyeun_norm, sample, any_of(RN2017_raw$EnsemblID)) # needed to not get errors when projecting (Ensembl version?)

### add Espinet IFNsign enrichment
rownames(Sauyeun_norm_) <- translate[rownames(Sauyeun_norm_)]
sample_info <- gsva(Sauyeun_norm_, IFNsign_geneset) %>% 
  t() %>% 
  as_tibble(rownames = "sample") %>% 
  right_join(sample_info, "sample")


## CPTAC
### processing
if (filter == TRUE){
  CPTAC_tumor <- dataFilter(CPTAC_tumor_raw, clinical_data, Molecular_phenotype_data, estimation_method, neoplastic_min, acinar_max, islet_max, tumor_tissue_min)
} else {
  CPTAC_tumor <- dplyr::select(CPTAC_tumor_raw, -low_purity_samples$case_id)
}

### Create sample_info_CPTAC

type_pamg <- projectMolGrad(newexp = column_to_rownames(dplyr::select(CPTAC_tumor,-EnsemblID), "Gene"),  geneSymbols = CPTAC_tumor$Gene) %>%
  as_tibble(rownames = "sample")

type_IFNsign <- gsva(data.matrix(column_to_rownames(dplyr::select(CPTAC_tumor,-EnsemblID), "Gene")), IFNsign_geneset) %>%
  t() %>% 
  as_tibble(rownames = "sample")

sample_info_CPTAC <- dplyr::select(type_pamg, sample, ICGCrnaseq) %>% 
  left_join(type_IFNsign, "sample") %>%
  dplyr::rename(PAMG = ICGCrnaseq)


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

type_IFNsign <- gsva(PACAOMICS_norm_, IFNsign_geneset) %>%
  t() %>% 
  as_tibble(rownames = "sample")

sample_info_PACAOMICS <- dplyr::select(sample_info, -PAMG, -IFNsign) %>%
  right_join(dplyr::select(type_pamg, sample, PDX)) %>% 
  left_join(type_IFNsign, "sample") %>%
  dplyr::rename(PAMG = PDX)


## CCLE
### Normalize (all together)
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

ccle_translate = deframe(ccle_raw[c("EnsemblID", "GeneName")])

ccle_norm_gn <- as_tibble(ccle_norm_, rownames = "EnsemblID") %>%
  mutate(GeneName = ccle_translate[EnsemblID]) %>% 
  dplyr::select(-EnsemblID) %>% group_by(GeneName) %>% 
  summarise_all(sum)


# PCA
## Sauyeun PDX sample and gene filtering
### get 5 and 5 most extreme samples for ICA3 
ICA3_density <- ggdensity(
  sample_info, x = "ICA3",
  rug = TRUE
)
ICA3_density

top_samples <- sample_info %>%
  arrange(ICA3) %>%
  dplyr::slice( unique(c(1:5, n() - 0:4)) ) %>%
  mutate(ISRact = ifelse(ICA3 < 0, "low_ICA3", "high_ICA3")) %>%
  arrange(sample) 

### get the most correlated genes with ICA3

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


top_genes <- mostCorrGenes(Sauyeun_norm, select_method, keepgenes, cormethod, top_samples, sample_info)

#### Check if they are the same as ISRact best subset 
table(top_genes$EnsemblID %in% colnames(ISRact_data_bestsubset)) #v different 568 genes missing

#### Check presence in datasets
table(top_genes$EnsemblID %in% CPTAC_tumor$EnsemblID) # of these 321 do not exist in the original

#### Final subset
Sauyeun_norm_subset <- dplyr::select(Sauyeun_norm, sample, top_genes$EnsemblID) %>% #to keep Jacobo's genes: any_of(colnames(ISRact_data_bestsubset))
  dplyr::filter(sample %in% top_samples$sample) %>%
  arrange(sample) 


## Run the PCA 
pca_pdx <- Sauyeun_norm_subset %>% 
  column_to_rownames( "sample") %>% 
  dplyr::select(any_of(CPTAC_tumor$EnsemblID)) %>%
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


## Projecting datasets on this PCA
### CPTAC
#### transpose human data for projection
human_data <- CPTAC_tumor %>%
  inner_join(dplyr::select(top_genes, EnsemblID), by = "EnsemblID") %>% #subset top genes
  pivot_longer(cols = 2:(length(CPTAC_tumor)-1), names_to = "case_id", values_to = "Expression") %>% 
  dplyr::filter(!str_detect(EnsemblID, 'NA')) %>% #remove eventual NA for gene Ensembl
  dplyr::select(-Gene) %>% 
  pivot_wider(names_from = EnsemblID, values_from = Expression, names_repair = "minimal") %>%
  column_to_rownames("case_id") 

# Project test: manual vs project()
projection_CPTAC <- scale(human_data[rownames(pca_pdx$rotation)], pca_pdx$center, pca_pdx$scale) %*% pca_pdx$rotation %>%
  as.data.frame()


p <- predict(pca_pdx, human_data) 
all(projection_CPTAC== p)


### PACAOMICS
projection_PACAOMICS <- predict(pca_pdx, PACAOMICS_norm) %>%
  as_tibble() %>%
  mutate(sample = PACAOMICS_norm$sample, .before = 1) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace(ISRact, is.na(ISRact) & sample %in% Sauyeun_norm$sample, "medium_ICA3"),
         ISRact = replace_na(ISRact, "Unknown"))

PACAOMICS_PC1 <- arrange(projection_PACAOMICS, PC1) %>% 
  dplyr::filter(!sample=="!") %>%
  mutate(PC1status = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3/3))), labels = c("low_PC1", "medium_PC1", "high_PC1"), include.lowest = TRUE)) %>%
  dplyr::select(sample, PC1,  PC1status) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace_na(ISRact, "medium_ICA3")) %>%
  inner_join(sample_info_PACAOMICS[, c("sample","PAMG", "ICA3")], by = "sample")


### Same original Sauyeun PDX with the middle samples
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
  inner_join(sample_info[, c("sample","PAMG", "ICA3", "IFNsign", "Diabetes")], by = "sample")


### CCLE
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
  left_join(ccle_info[,c("ccle_name","primary_tissue", "ISRact")], by="ccle_name")

ccle_PC1 <- arrange(projection_ccle, PC1) %>% 
  mutate(PC1status = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3/3))), labels = c("low_PC1", "medium_PC1", "high_PC1"), include.lowest = TRUE)) %>%
  dplyr::select(ccle_name, PC1,  PC1status) %>%
 


#Plot pca and projections
#Add ISR status to PCA df
pca_full_df <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  inner_join(top_samples[,c("sample","ISRact")])

#-------------------------------------------------------------------------------
# Plot projecion PacaOmics and then projection Proteogenomics

show_projection <- bind_rows(as_tibble(pca_full_df) %>% mutate(dataset = "Sauyeun PDX"),
                             mutate(projection_PACAOMICS, sample = paste0(sample, "_OG"),
                                    dataset = "PACAOMICS PDX"),
                             as_tibble(projection_CPTAC) %>% mutate(ISRact = "Unknown",
                                                              dataset = "CPTAC"),
                             as_tibble(projection_ccle %>% mutate(ISRact = if_else(is.na(ISRact), "Unknown", ISRact),
                                                                  dataset = "CCLE")))

dplyr::mutate(show_projection, ISRact = str_replace(ISRact, 'ICA3', 'ISRact')) %>%
  dplyr::filter(dataset %in% c("Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
              ISRact %in% c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) %>% 
ggplot(aes(x=PC1, y=PC2, color = ISRact, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Sauyeun PDX", "PACAOMICS PDX", "CPTAC", "CCLE")) +
  scale_color_discrete(limits = c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) + #c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue))
  ylim(-250,15)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# CCLE plot with names and colored according to any gene expression
pivot_longer(ccle_norm_gn,-GeneName, names_to = "ccle_name", values_to = "expression") %>%
  pivot_wider(id_cols = "ccle_name", names_from = "GeneName", values_from = "expression") %>%
  inner_join(projection_ccle, by="ccle_name") %>% 
  ggplot(aes(x=PC1, y=GATA6, label = str_remove(ccle_name, "_PANCREAS"))) +
  geom_point(aes(color = ISRact)) + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Distribution of PC1 gene weights with respect to IFNsign
PC1_vs_IFNsign <- as_tibble(pca_pdx$rotation, rownames = "EnsemblID") %>%
  mutate(Gene_symbol = translate[EnsemblID],
         in_IFNsign = ifelse(Gene_symbol %in% IFNsign_geneset$IFNsign, TRUE, FALSE))

ggplot(PC1_vs_IFNsign, aes(x = PC1)) +
  geom_density() +
  geom_rug(data = dplyr::filter(PC1_vs_IFNsign, in_IFNsign == T))

dplyr::filter(PC1_vs_IFNsign, in_IFNsign == T) %>% dplyr::select(PC1, Gene_symbol)
#-------------------------------------------------------------------------------
#Select human cohort most extreme samples regarding PC1 #! Remove when further scripts are adapted!
human_ISR <- arrange(projection_CPTAC, PC1) %>%
  dplyr::slice( unique(c(1:(nrow(projection_CPTAC)/3), n() - 0:((nrow(projection_CPTAC)/3)-1))) ) %>%
  mutate(ISRact = ifelse(PC1 < median(projection_CPTAC$PC1), "low_ISR", "high_ISR")) %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, ISRact)
#write.csv(human_ISR, "human_ISR1.csv", row.names = FALSE) #! Old! from Enzo's Code

# Select most extreme samples for further comparisons
# PACAOMICS
PACAOMICS_PC1 <-  arrange(projection_PACAOMICS, PC1) %>%
  mutate(PC1status = if_else(PC1 < quantile(projection_PACAOMICS$PC1, probs = 0.3333), "low_PC1",
                             if_else(PC1 < quantile(projection_PACAOMICS$PC1, probs = 0.6666), "medium_PC1", "high_PC1"))) %>% 
  dplyr::select(sample, PC1,  PC1status) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace_na(ISRact, "medium_ICA3")) %>%
  left_join(sample_info_PACAOMICS[, c("sample","PAMG", "IFNsign", "ICA3")], by = "sample")

write.csv(PACAOMICS_PC1, "../02_Output/PACAOMICS_PC1.csv", row.names = FALSE)

# CPTAC
CPTAC_PC1 <- arrange(projection_CPTAC, PC1) %>%
  mutate(PC1status = if_else(PC1 < quantile(projection_CPTAC$PC1, probs = 0.3333), "low_PC1",
                             if_else(PC1 < quantile(projection_CPTAC$PC1, probs = 0.6666), "medium_PC1", "high_PC1"))) %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, PC1,  PC1status) %>%
  left_join(sample_info_CPTAC, by = "sample")



write.csv(CPTAC_PC1, "../02_Output/CPTAC_PC1.csv", row.names = FALSE)


# Exploration of projection of CPTAC
## General Boxplot CPTAC clinical data
projection_full_df <- rownames_to_column(projection_CPTAC, "case_id") %>%
  inner_join(clinical_data[,c("case_id","medical_condition","tumor_stage_pathological")], by = "case_id" ) %>%
  inner_join(dplyr::select(Molecular_phenotype_data, case_id, Bailey, Collisson, Moffitt, contains("deconv"), contains("histology"), contains("pathway"), contains("MCPCounter")), by = "case_id")

projection_cor <- dplyr::select(projection_full_df, PC1:PC10, immune_deconv:`Myeloid dendritic cells_MCPCounter`) %>% #select continuous variables
  as.matrix() %>%
  Hmisc::rcorr() %>%
  map( ~data.frame(.x))


#function to format correlations in order to plot with ggplot2
formatted_cors <- function(df){
  df %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    dplyr::rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
}

formatted_cors(projection_cor) %>%
  dplyr::filter(str_detect(measure1, 'PC1\\b') & !str_detect(measure2, 'PC\\d')) %>%
  mutate(measure1 = fct_relevel(measure1, unique(measure1))) %>%
  mutate(measure2 = fct_relevel(measure2, unique(measure2))) %>%
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title="Correlations in projection",
       subtitle="Only significant Pearson's correlation coefficients shown") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  theme(text=element_text(family="Roboto"),
        axis.text.x = element_text(angle = 90, hjust=1))


## Anova between PC1 and cofactors (Diabete and Moffit)
projection_full_df$Diabetes <- ifelse(str_detect(projection_full_df$medical_condition, "Diabetes"), "yes", "no")

#Chose Moffitt or Diabetes for model
model  <- lm(PC1 ~ Moffitt, data = projection_full_df)
#model  <- lm(PC1 ~ Diabetes, data = projection_full_df)

#Test normality of residuals
ggqqplot(residuals(model))
shapiro_test(residuals(model))

#Test homogeneity of variances
plot(model, 1)
projection_full_df %>% levene_test(PC1 ~ Moffitt)

#Perform anova
res.aov <- projection_full_df %>% anova_test(PC1 ~ Moffitt)
res.aov

# Comparisons by paires
pwc <- projection_full_df %>% tukey_hsd(PC1 ~ Moffitt)
pwc

# Visualization : Boxplots with p-values
pwc <- pwc %>% add_xy_position(x = "Moffitt")
ggboxplot(projection_full_df, x = "Moffitt", y = "PC1") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )
#-------------------------------------------------------------------------------
# Custom better looking dotplot
as_tibble(projection_full_df) %>% ggplot(aes(y = PC1, x = Diabetes)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, alpha = 0.5) +
  geom_boxplot(aes(fill = Diabetes), position=position_dodge(0.8), width=0.05) +
  theme_classic()

as_tibble(projection_full_df) %>% ggplot(aes(y = PC1, x = Moffitt)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1, alpha = 0.5) +
  geom_boxplot(aes(fill = Moffitt), position=position_dodge(0.8), width=0.05) +
  scale_fill_manual(values = c("#FF7F00","#0000FC","grey")) +
  theme_classic()
#-------------------------------------------------------------------------------

# PC1 is Basal/Classical vs PC1 is ISRact
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

## Sauyeun PDX
### ISR vs PC1
correlation_plotter(data = Sauyeun_PC1, col1 = "ICA3", col2 = "PC1", data_name = "Sauyeun PDX")

### PAMG vs PC1
correlation_plotter(data = Sauyeun_PC1, col1 = "PAMG", col2 = "PC1", data_name = "Sauyeun PDX")

### ISR vs PAMG
correlation_plotter(data = Sauyeun_PC1, col1 = "ICA3", col2 = "PAMG", data_name = "Sauyeun PDX")

### ISR vs IFNsign
correlation_plotter(data = Sauyeun_PC1, col1 = "ICA3", col2 = "IFNsign", data_name = "Sauyeun PDX")

### PC1 vs IFNsign
correlation_plotter(data = Sauyeun_PC1, col1 = "PC1", col2 = "IFNsign", data_name = "Sauyeun PDX")


## PACAOMICs
### ISR vs PC1
correlation_plotter(data = PACAOMICS_PC1, col1 = "ICA3", col2 = "PC1", data_name = "PACAOMICs PDX")

### PAMG vs PC1
correlation_plotter(data = PACAOMICS_PC1, col1 = "PAMG", col2 = "PC1", data_name = "PACAOMICs PDX")

### ISR vs PAMG
correlation_plotter(data = PACAOMICS_PC1, col1 = "ICA3", col2 = "PAMG", data_name = "PACAOMICs PDX")

### ISR vs IFNsign
correlation_plotter(data = PACAOMICS_PC1, col1 = "ICA3", col2 = "IFNsign", data_name = "PACAOMICs PDX")

### PC1 vs IFNsign
correlation_plotter(data = PACAOMICS_PC1, col1 = "PC1", col2 = "IFNsign", data_name = "PACAOMICs PDX")


## CPTAC
### PAMG vs PC1
correlation_plotter(data = CPTAC_PC1, col1 = "PAMG", col2 = "PC1", data_name = "CPTAC tumors")

### IFNsign vs PC1
correlation_plotter(data = CPTAC_PC1, col1 = "PC1", col2 = "IFNsign", data_name = "CPTAC tumors")


# Survival curves reguarding ISR status
surv_data <- dplyr::rename(human_ISR, case_id = sample) %>%
  inner_join(clinical_data, by = "case_id") 

fit <-   survfit(Surv(follow_up_days, status) ~ ISRact, data = surv_data)
print(fit)

# Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

