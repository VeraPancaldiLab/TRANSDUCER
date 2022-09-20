#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
source("functions.R")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/")

# PARAMETERS
#-------------------------------------------------------------------------------
component_reorientation = TRUE
reorient_cyt <- c(1, 1, 1, -1, 1, 1)
#-------------------------------------------------------------------------------

# load annotation with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "mmusculus_gene_ensembl",
                        version = 75)

#listAttributes(ensembl75, page="feature_page")
annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id',
                                        'entrezgene',
                                        'mgi_id',
                                        'chromosome_name'), mart = ensembl75)

gene_to_ensembl = deframe(annot_ensembl75[c( "external_gene_id", "ensembl_gene_id")])
ensembl_to_gene = setNames(names(gene_to_ensembl), gene_to_ensembl)

#Data loading
## ICs
ICA_cyt <- read_rds("02_Output/ICA_cyt.RDS")
ICA_pol <- read_rds("02_Output/ICA_pol.RDS")
ICA_TEs <- read_rds("02_Output/ICA_TEs.RDS")

if (component_reorientation == TRUE){
  ICA_cyt[["A"]] <- t(t(ICA_cyt[["A"]]) * reorient_cyt)
  ICA_cyt[["S"]] <- t(t(ICA_cyt[["S"]]) * reorient_cyt)
}

## Norm_data
cyt_m <- read_tsv("../../00_Data/Processed_data/normHost_Cyt.tsv")
pol_m <- read_tsv("../../00_Data/Processed_data/normHost_Pol.tsv")

## Metadata 
annotations <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>%
  dplyr::select(sample, PAMG, ICA1, ICA3) %>%
  dplyr::rename(ISRact = ICA3)
  
# Comparison of sample weights
## component wide
A_cyt <- as_tibble(ICA_cyt$A, rownames = "sample") %>% rename_with( ~paste0(.,".cyt"), -sample)
A_pol <- as_tibble(ICA_pol$A, rownames = "sample") %>% rename_with( ~paste0(.,".pol"), -sample)
A_TEs <- as_tibble(ICA_TEs$A, rownames = "sample") %>% rename_with( ~paste0(.,".TEs"), -sample)

Acomparison <- inner_join(A_cyt, A_pol, by = "sample") %>%
  inner_join(A_TEs, by = "sample") %>%
  inner_join(annotations, by = "sample") %>%
  column_to_rownames("sample")

### cyt ~ pol
dplyr::select(Acomparison, ends_with(c("cyt", "pol", "ISRact", "PAMG"))) %>%
  formatted_cors(cor.stat = "pearson") %>%
  dplyr::filter(str_detect(measure1, ".pol", negate = TRUE), str_detect(measure2, ".cyt", negate = TRUE)) %>%
  ggplot(aes(measure2, measure1, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "Cytosolic vs Polysome", 
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  rotate_x_text(angle = 45)

### cyt ~ TEs
dplyr::select(Acomparison, ends_with(c("cyt", "TEs", "ISRact", "PAMG"))) %>%
    formatted_cors(cor.stat = "pearson") %>%
  dplyr::filter(str_detect(measure1, ".TEs", negate = TRUE), str_detect(measure2, ".cyt", negate = TRUE)) %>%
  ggplot(aes(measure2, measure1, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "Cytosolic vs Translation Efficacies", 
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
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "Polysome vs Translation Efficacies", 
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
  labs(x = NULL, y = NULL, fill = "Pearson's\nCorrelation", title= "Cytosolic vs Polysome (gene weights)", 
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
column_annot <- dplyr::select(A_TEs, sample, IC.3.TEs, IC.4.TEs, IC.6.TEs) %>%
  inner_join(annotations, by = "sample") %>% column_to_rownames("sample")

order_by <- arrange(column_annot, ISRact) %>% rownames()
clust = F

## Translation initiation
translation_initiation <- dplyr::filter(reactome, gs_name == "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION")
ensembl_genesym <- dplyr::select(translation_initiation , c(ensembl_gene, gene_symbol)) %>% deframe()

dplyr::filter(cyt_m, EnsemblID %in% translation_initiation$ensembl_gene) %>%
  mutate(gene_symbol = ensembl_genesym[EnsemblID]) %>%
  dplyr::select(gene_symbol, A_cyt$sample) %>% 
  relocate(gene_symbol, all_of(order_by)) %>%
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

# CAF subtype markers

## Elayada et al 2019
### gene weight check
Elayada <- read_tsv("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/01_Input/CAF_Elayada_markers.tsv")
S_cyt_iCAFs <- dplyr::filter(S_cyt, EnsemblID %in% gene_to_ensembl[Elayada$iCAFs])
S_cyt_apCAFs <- dplyr::filter(S_cyt, EnsemblID %in% gene_to_ensembl[Elayada$apCAFs])
S_cyt_myCAFs <- dplyr::filter(S_cyt, EnsemblID %in% gene_to_ensembl[Elayada$myCAFs])



ggplot(S_cyt) +
  aes_string("IC.4.cyt") +
  geom_density() +
  geom_vline(xintercept = 0) +
  geom_rug(data = S_cyt_iCAFs, aes(IC.4.cyt, color = "iCAFs"), outside = T) +
  geom_rug(data = S_cyt_apCAFs, aes(IC.4.cyt, color = "apCAFs")) +
  geom_rug(data = S_cyt_myCAFs, aes(IC.4.cyt, color = "myCAFs"), sides = "t") +
  coord_cartesian(clip = 'off') +
  theme_classic()

### heatmap of gene expression
annot_row <- pivot_longer(Elayada,  cols = tidyselect::everything(), names_to = "signature", values_to = "genes") %>%
  dplyr::filter(!is.na(genes)) %>%
  dplyr::arrange(genes) %>%
  mutate(value = 1) %>% 
  pivot_wider(names_from = "signature", id_cols = "genes") %>%
  replace(is.na(.), 0) %>% column_to_rownames("genes")

annot_col <- inner_join(A_cyt, annotations, "sample") %>% 
  dplyr::select(sample, PAMG, ISRact, IC.4.cyt) %>% column_to_rownames("sample")

annot_colors <- list(PAMG = c("#FF7F00", "white", "#377DB8"),
                     ISRact = c("#FFFFCC", "#006837"),
                     IC.4.cyt = c("#FFFFCC", "#5b0066")) 

stopifnot(names(cyt_m)[-1] == annot_col$sample)
mutate(cyt_m, Genenames = ensembl_to_gene[cyt_m$EnsemblID]) %>%
  dplyr::filter(Genenames %in% flatten_chr(Elayada), !is.na(Genenames)) %>%
  dplyr::select(!EnsemblID) %>%
  relocate(Genenames) %>% column_to_rownames("Genenames") %>% 
  pheatmap(scale = "row", color = colorRampPalette(c("#0C6291", "#FBFEF9", "#A63446"))(100),
           annotation_row = annot_row,
           annotation_col = annot_col,
           annotation_colors = annot_colors,
           cluster_rows = T, show_rownames = TRUE,
           show_colnames= FALSE, main = "Stromal transcription of Elayada et al. 2019 murine CAF markers")

## Verginadis et al 2022
### gene weight check
Verginadis <- read_tsv("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/01_Input/CAF_Verginadis_markers.tsv")
S_cyt_CAFs <- dplyr::filter(S_cyt, EnsemblID %in% gene_to_ensembl[Verginadis$CAFs])
S_cyt_vCAFs <- dplyr::filter(S_cyt, EnsemblID %in% gene_to_ensembl[Verginadis$vCAFs])
S_cyt_mCAFs <- dplyr::filter(S_cyt, EnsemblID %in% gene_to_ensembl[Verginadis$mCAFs])
S_cyt_cCAFs <- dplyr::filter(S_cyt, EnsemblID %in% gene_to_ensembl[Verginadis$cCAFs])
S_cyt_melCAFs <- dplyr::filter(S_cyt, EnsemblID %in% gene_to_ensembl[Verginadis$melCAFs])


ggplot(S_cyt) +
  aes_string("IC.4.cyt") +
  geom_density() +
  geom_vline(xintercept = 0) +
  geom_rug(data = S_cyt_vCAFs, aes(IC.4.cyt, color = "vCAFs"), outside = T) +
  geom_rug(data = S_cyt_mCAFs, aes(IC.4.cyt, color = "mCAFs")) +
  geom_rug(data = S_cyt_cCAFs, aes(IC.4.cyt, color = "cCAFs")) +
  geom_rug(data = S_cyt_melCAFs, aes(IC.4.cyt, color = "melCAFs"), outside = T) +
  geom_rug(data = S_cyt_CAFs, aes(IC.4.cyt, color = "CAFs"), sides = "t") +
  coord_cartesian(clip = 'off') +
  theme_classic()

### heatmap of gene expression
annot_row <- pivot_longer(Verginadis,  cols = tidyselect::everything(), names_to = "signature", values_to = "genes") %>%
  dplyr::filter(!is.na(genes)) %>%
  dplyr::arrange(genes) %>%
  mutate(value = 1) %>% 
  pivot_wider(names_from = "signature", id_cols = "genes") %>%
  replace(is.na(.), 0) %>% column_to_rownames("genes")

annot_col <- inner_join(A_cyt, annotations, "sample") %>% 
  dplyr::select(sample, PAMG, ISRact, IC.4.cyt) %>% column_to_rownames("sample")

annot_colors <- list(PAMG = c("#FF7F00", "white", "#377DB8"),
                     ISRact = c("#FFFFCC", "#006837"),
                     IC.4.cyt = c("#FFFFCC", "#5b0066")) 

stopifnot(names(cyt_m)[-1] == annot_col$sample)
mutate(cyt_m, Genenames = ensembl_to_gene[cyt_m$EnsemblID]) %>%
  dplyr::filter(Genenames %in% flatten_chr(Verginadis), !is.na(Genenames)) %>%
  dplyr::select(!EnsemblID) %>%
  relocate(Genenames) %>% column_to_rownames("Genenames") %>% 
  pheatmap(scale = "row", color = colorRampPalette(c("#0C6291", "#FBFEF9", "#A63446"))(100),
           annotation_row = annot_row,
           annotation_col = annot_col,
           annotation_colors = annot_colors,
           cluster_rows = T, show_rownames = TRUE,
           show_colnames= False, main = "Stromal transcription of Verginadis et al. 2022 murine CAF markers")
