#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
source("functions.R")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/")

#Data loading

## ICs
ICA_cyt <- read_rds("02_Output/ICA_cyt.RDS")
ICA_pol <- read_rds("02_Output/ICA_pol.RDS")
ICA_TEs <- read_rds("02_Output/ICA_TEs.RDS")

## Norm_data
cyt_m <- read_tsv("../../00_Data/Processed_data/normHost_Cyt.tsv")
pol_m <- read_tsv("../../00_Data/Processed_data/normHost_Pol.tsv")

## Metadata 
cancer_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>%
  dplyr::select(sample, PAMG, ICA1, ICA3) %>%
  dplyr::rename(ISRact = ICA3)
  
# Comparison of sample weights
## component wide
A_cyt <- as_tibble(ICA_cyt$A, rownames = "sample") %>% rename_with( ~paste0(.,".cyt"), -sample)
A_pol <- as_tibble(ICA_pol$A, rownames = "sample") %>% rename_with( ~paste0(.,".pol"), -sample)
A_TEs <- as_tibble(ICA_TEs$A, rownames = "sample") %>% rename_with( ~paste0(.,".TEs"), -sample)

Acomparison <- inner_join(A_cyt, A_pol, by = "sample") %>%
  inner_join(A_TEs, by = "sample") %>%
  inner_join(cancer_info, by = "sample") %>%
  column_to_rownames("sample")

### cyt ~ pol
dplyr::select(Acomparison, ends_with(c("cyt", "pol"))) %>%
  formatted_cors(cor.stat = "pearson") %>%
  dplyr::filter(str_detect(measure1, ".cyt"), str_detect(measure2, ".pol")) %>%
  ggplot(aes(measure2, measure1, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "cyt vs pol", 
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  rotate_x_text(angle = 45)

### cyt ~ TEs
dplyr::select(Acomparison, ends_with(c("cyt", "TEs"))) %>%
    formatted_cors(cor.stat = "pearson") %>%
  dplyr::filter(str_detect(measure1, ".cyt"), str_detect(measure2, ".TEs")) %>%
  ggplot(aes(measure2, measure1, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "cyt vs TEs", 
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  rotate_x_text(angle = 45)

### pol ~ TEs
dplyr::select(Acomparison, ends_with(c("pol", "TEs"))) %>%
  formatted_cors(cor.stat = "pearson") %>%
  dplyr::filter(str_detect(measure1, ".pol"), str_detect(measure2, ".TEs")) %>%
  ggplot(aes(measure2, measure1, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "pol vs TEs", 
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  rotate_x_text(angle = 45)

## gene wise
S_cyt <- as_tibble(ICA_cyt$S, rownames = "EnsemblID") %>% rename_with( ~paste0(.,".cyt"), -EnsemblID)
S_pol <- as_tibble(ICA_pol$S, rownames = "EnsemblID") %>% rename_with( ~paste0(.,".pol"), -EnsemblID)
S_TEs <- as_tibble(ICA_TEs$S, rownames = "EnsemblID") %>% rename_with( ~paste0(.,".TEs"), -EnsemblID)

Scomparison <- inner_join(S_cyt, S_pol, by = "EnsemblID") %>%
  inner_join(S_TEs, by = "EnsemblID") %>%
  column_to_rownames("EnsemblID")

### cyt ~ pol
dplyr::select(Scomparison, ends_with(c("cyt", "pol"))) %>%
  formatted_cors(cor.stat = "pearson") %>%
  dplyr::filter(str_detect(measure1, ".cyt"), str_detect(measure2, ".pol")) %>%
  ggplot(aes(measure2, measure1, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "cyt vs pol (gene weights)", 
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  rotate_x_text(angle = 45)


# Reactome pathways cyt/pol vs. ICA.TEs
all_genesets <- msigdbr("Mus musculus")
reactome <- dplyr::filter(all_genesets, gs_subcat == "CP:REACTOME")
column_annot <- dplyr::select(A_TEs, sample, IC.3.TEs, IC.4.TEs, IC.7.TEs) %>%
  inner_join(cancer_info, by = "sample") %>% column_to_rownames("sample")

order_by <- arrange(column_annot, ISRact) %>% rownames()
clust = F

## Translation initiation
translation_initiation <- dplyr::filter(reactome, gs_name == "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION")
ensembl_genesym <- dplyr::select(translation_initiation , c(ensembl_gene, gene_symbol)) %>% deframe()

dplyr::filter(cyt_m, EnsemblID %in% translation_initiation$ensembl_gene) %>%
  mutate(gene_symbol = ensembl_genesym[EnsemblID]) %>%
  dplyr::select(gene_symbol, A_cyt$sample) %>% 
  relocate(gene_symbol, order_by) %>%
  column_to_rownames("gene_symbol") %>% 
  pheatmap(scale = "row", annotation_col = column_annot,
           cluster_cols = clust,
           cluster_rows = T, fontsize = 5, 
           treeheight_col = 10, main = "T.initiation cyt vs ICA.TEs", cellheight = 5, filename = "02_Output/reactome_tinitiation.png")

## Translation elongation
translation_elongation <- dplyr::filter(reactome, gs_name == "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION")
ensembl_genesym <- dplyr::select(translation_elongation , c(ensembl_gene, gene_symbol)) %>% deframe()

dplyr::filter(cyt_m, EnsemblID %in% translation_elongation$ensembl_gene) %>%
  mutate(gene_symbol = ensembl_genesym[EnsemblID]) %>%
  dplyr::select(gene_symbol, A_cyt$sample) %>% 
  relocate(gene_symbol, order_by) %>%
  column_to_rownames("gene_symbol") %>% 
  pheatmap(scale = "row", annotation_col = column_annot,
           cluster_cols = clust, 
           cluster_rows = T, fontsize = 5, 
           treeheight_col = 10, main = "T.elongation cyt vs ICA.TEs", cellheight = 5, filename = "02_Output/reactome_tenlongation.png")


## tRNA aminocydation
trna_aa <- dplyr::filter(reactome, gs_name == "REACTOME_TRNA_AMINOACYLATION")
ensembl_genesym <- dplyr::select(trna_aa , c(ensembl_gene, gene_symbol)) %>% deframe()

dplyr::filter(cyt_m, EnsemblID %in% trna_aa$ensembl_gene) %>%
  mutate(gene_symbol = ensembl_genesym[EnsemblID]) %>%
  dplyr::select(gene_symbol, A_cyt$sample) %>% 
  relocate(gene_symbol, order_by) %>%
  column_to_rownames("gene_symbol") %>% 
  pheatmap(scale = "row", annotation_col = column_annot,
           cluster_cols = clust, 
           cluster_rows = T, fontsize = 5, 
           treeheight_col = 10, main = "tRNA Aminoacylation cyt vs ICA.TEs", cellheight = 5, filename = "02_Output/reactome_aminoacylation.png")

## mTOR signaling
mtor <- dplyr::filter(reactome, gs_name == "REACTOME_MTOR_SIGNALLING")
ensembl_genesym <- dplyr::select(mtor , c(ensembl_gene, gene_symbol)) %>% deframe()

dplyr::filter(cyt_m, EnsemblID %in% mtor$ensembl_gene) %>%
  mutate(gene_symbol = ensembl_genesym[EnsemblID]) %>%
  dplyr::select(gene_symbol, A_cyt$sample) %>% 
  relocate(gene_symbol, order_by) %>%
  column_to_rownames("gene_symbol") %>% 
  pheatmap(scale = "row", annotation_col = column_annot,
           cluster_cols = clust, 
           cluster_rows = T, fontsize = 5, 
           treeheight_col = 10, main = "MTOR signaling cyt vs ICA.TEs", cellheight = 5, filename = "02_Output/reactome_mTOR.png")