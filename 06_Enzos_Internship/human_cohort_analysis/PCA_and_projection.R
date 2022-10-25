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
mRNA_tumor <- read_delim("PDAC_LinkedOmics_Data/mRNA_RSEM_UQ_log2_Tumor.cct", 
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
Human_Tumor_rawcount_Transcriptome <- read_delim("Human-Tumor_rawcount_Transcriptome.tsv", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)

### Extended 73 sample cohort


################################################################################
# PARAMETERS
rawdata = rawTumor_Cyt #geneCount_raw_28s_totalRNA | rawTumor_Cyt
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

#Add a Gene names column to rawdata
translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

rawdata$EnsemblID %>% 
  translate[.] %>%
  make.names(unique = TRUE) -> rawdata$Genenames


#Add a Ensemble ID column to mRNA_tumor (Ensembl release 95 according to the paper)
#https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=genes&hgta_track=mane&hgta_table=mane&hgta_doSchema=describe+table+schema

reverse_translate = deframe(annot_ensembl95[c("external_gene_name", "ensembl_gene_id")])

mRNA_tumor$Gene %>% 
  reverse_translate[.] %>%
  make.names(unique = TRUE) -> mRNA_tumor$EnsemblID


# Processing and Normalization
## Sauyeun PDX
### Normalize raw counts
normdata <- dplyr::select(rawdata, -Genenames) %>%
  column_to_rownames( "EnsemblID") %>%
  DGEList() %>%
  calcNormFactors(method= norm_method) %>%
  cpm(log=TRUE) %>%
  t() %>%
  as_tibble(rownames = "sample")

### Visualize it
data <- pivot_longer(normdata, cols = 2:length(normdata), names_to = "EnsemblID", values_to = "Expression", names_repair = "unique")

g1 <- ggplot(data , mapping = aes(x = sample, y = Expression, fill = sample))+
  geom_violin()+
  theme(legend.position = "none")+
  coord_flip()

g1

### Keep common genes with Human_Tumor_rawcount_Transcriptome 
normdata <- dplyr::select(normdata, sample, any_of(Human_Tumor_rawcount_Transcriptome$EnsemblID)) # needed to not get errors when projecting (Ensembl version?)



## CPTAC
### processing
if (filter == TRUE){
  mRNA_tumor <- dataFilter(mRNA_tumor, clinical_data, Molecular_phenotype_data, estimation_method, neoplastic_min, acinar_max, islet_max, tumor_tissue_min)
} else {
  mRNA_tumor <- dplyr::select(mRNA_tumor, -low_purity_samples$case_id)
}


## PACAOMICS
Human_Tumor_rawcount <- column_to_rownames(Human_Tumor_rawcount_Transcriptome, "EnsemblID") 

Human_Tumor_norm <- DGEList(Human_Tumor_rawcount) %>%
  calcNormFactors(method= norm_method) %>%
  cpm(log=TRUE) %>%
  t() %>%
  as_tibble(rownames = "sample")


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


top_genes <- mostCorrGenes(normdata, select_method, keepgenes, cormethod, top_samples, sample_info)

#### Check if they are the same as ISRact best subset 
table(top_genes$EnsemblID %in% colnames(ISRact_data_bestsubset)) #v different 568 genes missing

#### Check presence in datasets
table(top_genes$EnsemblID %in% mRNA_tumor$EnsemblID) # of these 321 do not exist in the original

#### Final subset
normdata_topg_s <- dplyr::select(normdata, sample, top_genes$EnsemblID) %>% #to keep Jacobo's genes: any_of(colnames(ISRact_data_bestsubset))
  dplyr::filter(sample %in% top_samples$sample) %>%
  arrange(sample) 


## Run the PCA 
pca_pdx <- normdata_topg_s %>% 
  column_to_rownames( "sample") %>% 
  dplyr::select(any_of(mRNA_tumor$EnsemblID)) %>%
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
human_data <- mRNA_tumor %>%
  inner_join(dplyr::select(top_genes, EnsemblID), by = "EnsemblID") %>% #subset top genes
  pivot_longer(cols = 2:(length(mRNA_tumor)-1), names_to = "case_id", values_to = "Expression") %>% 
  dplyr::filter(!str_detect(EnsemblID, 'NA')) %>% #remove eventual NA for gene Ensembl
  dplyr::select(-Gene) %>% 
  pivot_wider(names_from = EnsemblID, values_from = Expression, names_repair = "minimal") %>%
  column_to_rownames("case_id") 

# Project test: manual vs project()
projection <- scale(human_data[rownames(pca_pdx$rotation)], pca_pdx$center, pca_pdx$scale) %*% pca_pdx$rotation %>%
  as.data.frame()


p <- predict(pca_pdx, human_data) 
all(projection == p)


### PACAOMICS
projection_PACAOMICS <- predict(pca_pdx, Human_Tumor_norm) %>%
  as_tibble() %>%
  mutate(sample = Human_Tumor_norm$sample, .before = 1) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace_na(ISRact, "medium_ICA3"))

PACAOMICS_PC1 <- arrange(projection_PACAOMICS, PC1) %>% 
  dplyr::filter(!sample=="!") %>%
  mutate(PC1status = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3/3))), labels = c("low_PC1", "medium_PC1", "high_PC1"), include.lowest = TRUE)) %>%
  dplyr::select(sample, PC1,  PC1status) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace_na(ISRact, "medium_ICA3")) %>%
  inner_join(sample_info[, c("sample","PAMG", "ICA3")], by = "sample")


### Same original Sauyeun PDX with the middle samples
projection_Sauyeun <- predict(pca_pdx, normdata) %>%
  as_tibble() %>%
  mutate(sample = normdata$sample, .before = 1) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace_na(ISRact, "medium_ICA3")) 

Sauyeun_PC1 <- arrange(projection_Sauyeun, PC1) %>% 
  dplyr::filter(!sample=="!") %>%
  mutate(PC1status = cut(.$PC1, breaks = c(quantile(.$PC1, c(0:3/3))), labels = c("low_PC1", "medium_PC1", "high_PC1"), include.lowest = TRUE)) %>%
  dplyr::select(sample, PC1,  PC1status) %>%
  left_join(top_samples[,c("sample","ISRact")]) %>%
  mutate(ISRact = replace_na(ISRact, "medium_ICA3")) %>%
  inner_join(sample_info[, c("sample","PAMG", "ICA3", "Diabetes")], by = "sample")


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
                             as_tibble(projection) %>% mutate(ISRact = "Unknown",
                                                              dataset = "CPTAC"))

dplyr::mutate(show_projection, ISRact = str_replace(ISRact, 'ICA3', 'ISRact')) %>%
  dplyr::filter(dataset %in% c("Sauyeun PDX", "PACAOMICS PDX"),
              ISRact %in% c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) %>% 
ggplot(aes(x=PC1, y=PC2, color = ISRact, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Sauyeun PDX", "PACAOMICS PDX", "CPTAC")) +
  scale_color_discrete(limits = c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) +
  ylim(-200,15)
#-------------------------------------------------------------------------------

#Select human cohort most extreme samples regarding PC1  
human_ISR <- arrange(projection, PC1) %>%
  dplyr::slice( unique(c(1:(nrow(projection)/3), n() - 0:((nrow(projection)/3)-1))) ) %>%
  mutate(ISRact = ifelse(PC1 < median(projection$PC1), "low_ISR", "high_ISR")) %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, ISRact)

#Split in 3 group regarding projection dim 1
human_PC1 <- arrange(projection, PC1) %>%
  mutate(PC1status = cut(.$PC1, 3, labels = c("low_PC1", "medium_PC1", "high_PC1"))) %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, PC1,  PC1status)

#Write it on a csv to use it on another script
#write.csv(human_ISR, "human_ISR1.csv", row.names = FALSE)


# Exploration of projection of CPTAC
## General Boxplot CPTAC clinical data
projection_full_df <- rownames_to_column(projection, "case_id") %>%
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
## PacaOmics
paca_ICA3_PC1_cor <- cor(PACAOMICS_PC1$PC1, PACAOMICS_PC1$ICA3, method = "spearman")

ggplot(PACAOMICS_PC1, aes(PC1, ICA3, color = ICA3)) +
  geom_point(shape = 16, size = 4, show.legend = FALSE) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  labs(title = "Comparison between PC1 and ICA3 in PaCaOmics samples",
       subtitle = paste("Spearman correlation = ", round(paca_ICA3_PC1_cor, 2), sep = ""))

paca_PAMG_PC1_cor <- cor(PACAOMICS_PC1$PC1, PACAOMICS_PC1$PAMG, method = "spearman")

ggplot(PACAOMICS_PC1, aes(PC1, PAMG, color = PAMG)) +
  geom_point(shape = 16, size = 4, show.legend = FALSE) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e")  +
  labs(title = "Comparison between PC1 and PAMG in PaCaOmics samples",
       subtitle = paste("Spearman correlation = ", round(paca_PAMG_PC1_cor, 2), sep = ""))


## Sauyeun PDX
sauy_ICA3_PC1_cor <- cor(Sauyeun_PC1$PC1, Sauyeun_PC1$ICA3, method = "spearman")

ggplot(Sauyeun_PC1, aes(PC1, ICA3, color = ICA3)) +
  geom_point(shape = 16, size = 4, show.legend = FALSE) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  labs(title = "Comparison between PC1 and ICA3 in Sauyeun samples",
       subtitle = paste("Spearman correlation = ", round(sauy_ICA3_PC1_cor, 2), sep = ""))

sauy_PAMG_PC1_cor <- cor(Sauyeun_PC1$PC1, Sauyeun_PC1$PAMG, method = "spearman")

ggplot(Sauyeun_PC1, aes(PC1, PAMG, color = PAMG)) +
  geom_point(shape = 16, size = 4, show.legend = FALSE) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e")+
  labs(title = "Comparison between PC1 and PAMG in Sauyeun samples",
       subtitle = paste("Spearman correlation = ", round(sauy_PAMG_PC1_cor, 2), sep = ""))


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

