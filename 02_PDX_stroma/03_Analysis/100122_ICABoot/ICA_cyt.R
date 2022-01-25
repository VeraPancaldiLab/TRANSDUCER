#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
library(JADE)
library(corrplot)
library(Hmisc)
library(pheatmap)
library("ggpubr")
library("ggplotify")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/")
source("functions.R")
run_boot <- FALSE
range.comp <- 2:15 # when ncomp is =< df Warning: In sqrt(puiss[rangeW]) : NaNs produced
boot.iter <- 500

# load cyt data
cyt <- read_tsv("../../00_Data/Processed_data/normHost_Cyt.tsv")
sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>%
  column_to_rownames("sample") %>% .[colnames(cyt)[-1],]

sample_info$Diabetes <- as_factor(sample_info$Diabetes)
# load annotation with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "mmusculus_gene_ensembl",
                        version = 75)

#listAttributes(ensembl75, page="feature_page")
annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                          'external_gene_id',
                          'entrezgene',
                          'chromosome_name'), mart = ensembl75)

# filtering
## XY
annot_ensembl75 %>% dplyr::filter(!chromosome_name %in% c("X", "Y")) %>%
  pull(ensembl_gene_id) -> non_sex

cyt %>% dplyr::filter(EnsemblID %in% non_sex) -> cyt_

## mean centring gene wise
cyt_ %>% column_to_rownames("EnsemblID") %>%
  apply(1, function(x) x - mean(x)) %>% t() %>%
  data.frame() -> cyt__

## Inter quartile range (measure variability)
cyt__ %>% apply(1, IQR) -> iqrs
mostvar <- iqrs[iqrs > median(iqrs)]
cyt__ %>% rownames_to_column("EnsemblID") %>%
  dplyr::filter(EnsemblID %in% names(mostvar)) %>%
  column_to_rownames("EnsemblID") -> cyt_icaready

# celltype deconvolution (unfiltered data)


# ICA
## Bootstrapping
if (run_boot == TRUE){
  ### Baseline before bootstrap
  cyt_icaready %>% jade_range(range.comp, MARGIN = 1) -> base_res_gene
  cyt_icaready %>% jade_range(range.comp, MARGIN = 2) -> base_res_sample
  
  ### Bootstrap
  gene_boot <- jade_choosencom(cyt_icaready, base_res_gene,
                               MARGIN = 1,
                               iterations = boot.iter,
                               seed = 0
  )
  
  sample_boot <- jade_choosencom(cyt_icaready, base_res_sample,
                                 MARGIN = 2,
                                 iterations = boot.iter,
                                 seed = 0
  )
  
  boot_plots(s_boot = sample_boot, g_boot = gene_boot, name = "Cyt")
} 

## Most robust n.comp analysis
elected_ncomp <- 7 
#elected_ncomp <- 710
### Run best n.comp ICA
jade_result <- JADE(cyt_icaready, n.comp = elected_ncomp)
colnames(jade_result[["A"]]) <- paste("IC", 1:elected_ncomp, sep = ".")
rownames(jade_result[["A"]]) <- names(jade_result$Xmu)

### Component distribution analysis
A_mat <- as.data.frame(jade_result[["A"]])
S_mat <- as.data.frame(jade_result[["S"]])
annotations <- sample_info[-1]
stopifnot(rownames(A_mat) == rownames(annotations))
plot_sample_weights(A_mat, annotations)

### metadata
#### Corplot
corr_continuous <- annotations %>% dplyr::select(!Diabetes) %>% bind_cols(A_mat)
corr_continuous <- corr_continuous[rownames(A_mat),] # merge mess with the order


continuous_rcorr <- rcorr(data.matrix(corr_continuous), type = "spearman")
continuous_rcorr$r <- continuous_rcorr$r[colnames(annotations %>% dplyr::select(!Diabetes)),colnames(A_mat)]
continuous_rcorr$P <- continuous_rcorr$P[colnames(annotations %>% dplyr::select(!Diabetes)),colnames(A_mat)]

corrplot(continuous_rcorr$r,
         p.mat = continuous_rcorr$P, sig.level = 0.05, insig = "blank")
stopifnot(rownames(A_mat) == rownames(annotations))
bind_cols(A_mat, annotations[,c("ICA3", "PAMG")]) %>%
  rename(SerDep = ICA3) -> complete_annotation

### Supervised celltype deconvolution
mMCPcounter_res <-  read_tsv("../180122_Various/02_Output/mMCPcounter_results.tsv") %>% column_to_rownames("cell_types") %>% t() %>% as_tibble(rownames = "samples")
ImmuCC_res <-  read_tsv("../180122_Various/02_Output/ImmuCC_results.tsv")

# deconv <- mMCPcounter_res %>% column_to_rownames("samples") %>% .[rownames(A_mat),]
# decon_title <- "mMCPcounter"

deconv <- ImmuCC_res %>% column_to_rownames("samples") %>% .[rownames(A_mat),]
decon_title <- "ImmuCC"


#### Heatmap
##### Full (with clustering)
deconv %>% t() %>% pheatmap(main=decon_title, scale = "row", annotation_col = complete_annotation)

#### Corplot
stopifnot(rownames(complete_annotation)==rownames(deconv))
corr_decon <- deconv %>% bind_cols(complete_annotation)
corr_decon <- corr_decon[rownames(complete_annotation),] # merge mess with the order

corr_decon <- rcorr(data.matrix(corr_decon), type = "spearman")
corr_decon$r <- corr_decon$r[colnames(deconv),colnames(complete_annotation)]
corr_decon$P <- corr_decon$P[colnames(deconv),colnames(complete_annotation)]

corrplot(corr_decon$r,
         p.mat = corr_decon$P, sig.level = 0.05, insig = "blank")


### TF activity
tumour_tf <- read_tsv("../180122_Various/02_Output/TFact_tumor_PanCan.tsv")
stroma_tf <- read_tsv("../180122_Various/02_Output/TFact_stroma_Gtex.tsv")

tf_activity <- tumour_tf %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))
tf_title <- "tumour"

# tf_activity <- stroma_tf %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))
# tf_title <- "stroma"

mostvar_TF <- Get_mostvar(tf_activity, 50)
#### Heatmap
##### Full (with
tf_activity %>%
  pheatmap(main=paste( "TF activity", tf_title), scale = "row", show_rownames = FALSE, annotation_col = complete_annotation)

##### 50 most var
mostvar_TF %>%
  pheatmap(main=tf_title, scale = "row", annotation_col = complete_annotation)

#### 50 best TFs by component
PlotBestCorr(complete_annotation, tf_activity, 10, analysis_name = paste(tf_title, "_best_TFs_cyt", sep = ""))


### Head and Tail genes
#### translate to gene names
S_mat %>% as_tibble(rownames = "EnsemblIDs") -> S_gn_
translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])
S_mat %>% rownames() %>% translate[.] -> S_gn_$Genenames
S_gn_ %>% dplyr::select(!EnsemblIDs) %>%
  relocate(Genenames) -> S_gn

#### load normalized expression data for representation
genes_toplot_ <- read_tsv("../180122_Various/01_Input/HostCyt_foranalysis.tsv")
genes_toplot_$EnsemblIDs %>% translate[.] -> genes_toplot_$Genenames
genes_toplot_ %>% dplyr::select(!EnsemblIDs) %>%
  relocate(Genenames) -> genes_toplot

#### plot
PlotGeneWeights(S_mat, genes_toplot, 25, complete_annotation, analysis_name = "gene_weights_cyt")

