#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
library(JADE)
library(ggpubr)
library(scico)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/")
source("functions.R")
run_boot <- FALSE
range.comp <- 2:15 # when ncomp is =< df Warning: In sqrt(puiss[rangeW]) : NaNs produced
boot.iter <- 500

# load cyt data
cyt <- read_tsv("../../00_Data/Processed_data/normHost_Cyt.tsv")
sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>%
  dplyr::filter(sample %in% colnames(cyt)[-1]) %>%
  column_to_rownames("sample")

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

## Most robust n.comp ICA
elected_ncomp <- 7 
#elected_ncomp <- 710
### ICA
jade_result <- JADE(cyt_icaready, n.comp = elected_ncomp)
colnames(jade_result[["A"]]) <- paste("IC", 1:elected_ncomp, sep = ".")
rownames(jade_result[["A"]]) <- names(jade_result$Xmu)

### Component distribution analysis
# Sample distribution in ICs
A_mat <- as.data.frame(jade_result[["A"]])
annotations <- sample_info[-1]
for (ann in colnames(annotations)[2]){
  print(typeof(annotations[[ann]]))
  rug_aes <- annotations[[ann]]
  rug_name <- ann
  comps_plots <- lapply(colnames(A_mat), function(ic){
    p <- 
      ggplot(A_mat) +
      aes_string(ic) +
      geom_density() + 
      geom_rug(aes(color = rug_aes)) +
      theme_classic() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank())
    
    if(is.integer(rug_aes)) {
      p <- p  +
        scale_color_discrete()
      
    } else if(is.double(rug_aes)){
      p <- p  +
        scico::scale_color_scico(palette = "hawai")
    } 
    
    # if (!(ic %in% c("IC.1", paste("IC", 1+trunc(elected_ncomp/2), sep = ".")))) { # this is to control wich cells should have a ylabel
    #   p <- p +
    #     theme(axis.title.y = element_blank())
    # }
    p
    
  })
  print(ggarrange(plotlist = comps_plots, common.legend = T))
}


