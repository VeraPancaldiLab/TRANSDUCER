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

sample_info_ICGC <- dplyr::select(type_pamg, sample, ICGCarray) %>% 
  dplyr::rename(PAMG = ICGCarray) %>% 
  left_join(rownames_to_column(ICGC_clinical, "sample"), "sample")

# Projecting datasets on this PCA
## ICGC
### transpose human data for projection
human_data <- ICGC_ensembl %>%
  pivot_longer(cols = -EnsemblID, names_to = "case_id", values_to = "Expression") %>% 
  dplyr::filter(!str_detect(EnsemblID, 'NA')) %>% #remove eventual NA for gene Ensembl
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

projection_ICGC <-  predict(pca_pdx, human_data) %>% 
  as_tibble(rownames = "sample")


ICGC_PC1 <- arrange(projection_ICGC, PC1) %>%
  mutate(PC1status = if_else(PC1 < quantile(projection_ICGC$PC1, probs = 0.3333), "low_PC1",
                             if_else(PC1 < quantile(projection_ICGC$PC1, probs = 0.6666), "medium_PC1", "high_PC1"))) %>%
  dplyr::select(sample, PC1,  PC1status) %>%
  left_join(sample_info_ICGC, by = "sample")

#Plot pca and projections
#Add ISR status to PCA df
pca_full_df <- pca_pdx[["x"]] %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  inner_join(top_samples[,c("sample","ISRact")])

#-------------------------------------------------------------------------------
# Plot projecion PacaOmics and then projection Proteogenomics

show_projection <- bind_rows(as_tibble(pca_full_df) %>% mutate(dataset = "Sauyeun PDX"),
                             as_tibble(projection_ICGC) %>% mutate(ISRact = "Unknown",
                                                                    dataset = "ICGC"))

dplyr::mutate(show_projection, ISRact = str_replace(ISRact, 'ICA3', 'ISRact')) %>%
  dplyr::filter(dataset %in% c("Sauyeun PDX","PACAOMICS PDX", "ICGC"), # "Sauyeun PDX","PACAOMICS PDX", "CPTAC", "CCLE"
                ISRact %in% c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) %>% 
  ggplot(aes(x=PC1, y=PC2, color = ISRact, shape = dataset)) +
  geom_point() +
  scale_shape_discrete(limits = c("Sauyeun PDX", "PACAOMICS PDX", "ICGC")) +
  scale_color_discrete(limits = c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown')) +  #c('low_ISRact', 'high_ISRact', 'medium_ISRact', 'Unknown'))unique(projection_ccle$primary_tissue))
  geom_rect(xmin = min(projection_ICGC$PC1),
            xmax = quantile(projection_ICGC$PC1, probs = 0.3333), 
            ymin = min(projection_ICGC$PC2), ymax = max(projection_ICGC$PC2), linewidth = 0, fill = "red", alpha =0.002) +
  
  geom_rect(xmin = quantile(projection_ICGC$PC1, probs = 0.6666),
            xmax = max(projection_ICGC$PC1),
            ymin = min(projection_ICGC$PC2), ymax = max(projection_ICGC$PC2), linewidth = 0, fill = "green", alpha =0.002)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Plot comparisons with Basal/Classical and ISRact
## PAMG vs PC1
correlation_plotter(data = ICGC_PC1, col1 = "PAMG", col2 = "PC1", data_name = "ICGC")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Survival curves regarding ISR status
surv_data <- dplyr::rename(ICGC_PC1, case_id = sample)
## OS
### Kaplan Meyer 33up vs 33down
fit <-  survfit(Surv(OS, OS_event) ~ PC1status, 
                data = dplyr::filter(surv_data, PC1status != "medium_PC1"))
print(fit)

#### Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("tomato3", "seagreen"))

### Cox Proportional hazzards model
cox.mod <- coxph(Surv(OS, OS_event) ~ PC1, 
                 data = as.data.frame(surv_data))
#### assumption checking
##### Linearity
plot(predict(cox.mod),
     residuals(cox.mod, type = "martingale"),
     xlab = "fitted", ylab = "Martingale residuals",
     main = "Residuals Plot", las = 1)
abline(h=0)
lines(smooth.spline(predict(cox.mod),
                    residuals(cox.mod, type = "martingale")),
      col="red")

##### Proportional Hazzards
par(mfrow = c(1,1))
plot(cox.zph(cox.mod)) # if failed, Proportional
abline(h=0, col =2)

#### Model plotting
ggforest(cox.mod)
ggadjustedcurves(cox.mod, data=as.data.frame(surv_data))


## PFS
### Kaplan Meyer 33up vs 33down
fit <-  survfit(Surv(PFS, PFS_event) ~ PC1status, 
                data = dplyr::filter(surv_data, PC1status != "medium_PC1"))
print(fit)

#### Change color, linetype by strata, risk.table color by strata
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("tomato3", "seagreen"))

### Cox Proportional hazzards model
cox.mod <- coxph(Surv(PFS, PFS_event) ~ PC1, 
                 data = as.data.frame(surv_data))
#### assumption checking
##### Linearity
plot(predict(cox.mod),
     residuals(cox.mod, type = "martingale"),
     xlab = "fitted", ylab = "Martingale residuals",
     main = "Residuals Plot", las = 1)
abline(h=0)
lines(smooth.spline(predict(cox.mod),
                    residuals(cox.mod, type = "martingale")),
      col="red")

##### Proportional Hazzards
par(mfrow = c(1,1))
plot(cox.zph(cox.mod)) # if failed, Proportional
abline(h=0, col =2)

#### Model plotting
ggforest(cox.mod)
ggadjustedcurves(cox.mod, data=as.data.frame(surv_data))
