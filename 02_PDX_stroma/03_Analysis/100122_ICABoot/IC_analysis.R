#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
library(ggVennDiagram)
library(msigdbr)
library(pheatmap)
library(Hmisc)
library(ggrepel)
source("functions.R")
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/")

# PARAMETERS
#-------------------------------------------------------------------------------
component_reorientation = TRUE
reorient_cyt <- c(1, 1, 1, -1, 1, 1)
reorient_TEs <- c(1, 1, -1, 1, 1, -1)
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
  ICA_TEs[["A"]] <- t(t(ICA_TEs[["A"]]) * reorient_TEs)
  ICA_TEs[["S"]] <- t(t(ICA_TEs[["S"]]) * reorient_TEs)
}

## Norm_data
cyt_m <- read_tsv("../../00_Data/Processed_data/normHost_Cyt.tsv")
pol_m <- read_tsv("../../00_Data/Processed_data/normHost_Pol.tsv")
TEs_m <- read_tsv("../081221_TranslationEfficacy/02_Output/TEs.tsv")

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

## 3d Scatterplot comparison
library(plotly)
axx <- list(
  title = "ISRact"
)

axy <- list(
  title = "IC.4 Total RNA"
)

axz <- list(
  title = "IC.6 TE"
)

colorscale <- list(c(0, 1), c("#FF7F00", "#377DB8"))

fig <- plot_ly(rownames_to_column(Acomparison, "sample"), x=~ISRact, y=~IC.4.cyt, z=~IC.6.TEs, text =~sample, 
               marker = list(color = ~PAMG, colorscale = colorscale, showscale = TRUE))
fig <- fig %>% add_markers() 
fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz),
                      annotations = list(x = 1.09,
                                         y = 1.02,
                                         text = 'Basal/Classical',
                                         xref = 'paper',
                                         yref = 'paper',
                                         showarrow = FALSE))
fig
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

# Reactome pathways genes vs enrichment
## preparation
all_genesets <- msigdbr("Mus musculus")
reactome <- dplyr::filter(all_genesets, gs_subcat == "CP:REACTOME")
annot_col <- inner_join(A_cyt, annotations, "sample") %>% 
  dplyr::select(sample, PAMG, ISRact, IC.4.cyt) %>% column_to_rownames("sample")

annot_colors <- list(PAMG = c("#FF7F00", "white", "#377DB8"),
                     ISRact = c("#FFFFCC", "#006837"),
                     IC.4.cyt = c("#FFFFCC", "#5b0066")) 

## plot code
PlotGenesetExpression <- function(expression, gene_set, title, filename, ...){
  dplyr::filter(expression, EnsemblID %in% gene_set$ensembl_gene) %>%
    mutate(gene_symbol = ensembl_to_gene[EnsemblID]) %>%
    dplyr::select(gene_symbol, A_cyt$sample) %>% 
    relocate(gene_symbol, all_of(order_by)) %>%
    column_to_rownames("gene_symbol") %>% 
    pheatmap(annotation_col = annot_col, 
             annotation_colors = annot_colors, fontsize = 5, 
             treeheight_col = 10, main = title, cellheight = 5, filename = filename, ...)
}

## mRNA activation by CAP and 43S binding
mRNAactivation <- dplyr::filter(reactome, str_detect(gs_name, "REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP*")) %>% 
  dplyr::select(gene_symbol, ensembl_gene)
title <- "mRNA activation by CAP and 43S binding cyt vs IC.4.cyt"
filename <- "02_Output/IC4cyt_mRNAactivation.png"
PlotGenesetExpression(cyt_m, mRNAactivation, title, filename)

## Translation initiation
translation_initiation <- dplyr::filter(reactome, gs_name == "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION") %>% 
  dplyr::select(gene_symbol, ensembl_gene)
title <- "T.initiation cyt vs IC.4.cyt"
filename <- "02_Output/IC4cyt_translation_initiation.png"
PlotGenesetExpression(cyt_m, translation_initiation, title, filename)

## Translation elongation
translation_elongation <- dplyr::filter(reactome, gs_name == "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION") %>% 
  dplyr::select(gene_symbol, ensembl_gene)
title <- "T.elongation cyt vs IC.4.cyt"
filename <- "02_Output/IC4cyt_translation_elongation.png"
PlotGenesetExpression(cyt_m, translation_elongation, title, filename)

## GCN2 Response to AA deprivation
GCN2AAdeprivation <- dplyr::filter(reactome, str_detect(gs_name, "REACTOME_RESPONSE_OF_EIF2AK4_GCN2")) %>% 
  dplyr::select(gene_symbol, ensembl_gene)

title <- "GCN2 response to AA deprivation cyt vs IC.4.cyt"
filename <- "02_Output/IC4cyt_GCN2AAdeprivation.png"
PlotGenesetExpression(cyt_m, GCN2AAdeprivation, title, filename)

## Translation
translation <- dplyr::filter(reactome, gs_name == "REACTOME_TRANSLATION") %>% 
  dplyr::select(gene_symbol, ensembl_gene)

title <- "Translation cyt vs IC.4.cyt"
filename <- "02_Output/IC4cyt_translation.png"
PlotGenesetExpression(cyt_m, translation, title, filename)

## SRP cotranslational targeting
SRPcotranstargeting <- dplyr::filter(reactome, gs_name == "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE") %>% 
  dplyr::select(gene_symbol, ensembl_gene)

title <- "SRP cotranslational targeting cyt vs IC.4.cyt"
filename <- "02_Output/IC4cyt_SRPcotranstargeting.png"
PlotGenesetExpression(cyt_m, SRPcotranstargeting, title, filename)

## rRNA processing 
rRNAprocessing <- dplyr::filter(reactome, gs_name == "REACTOME_RRNA_PROCESSING") %>% 
  dplyr::select(gene_symbol, ensembl_gene)

title <- "rRNA processing cyt vs IC.4.cyt"
filename <- "02_Output/IC4cyt_rRNAprocessing.png"
PlotGenesetExpression(cyt_m, rRNAprocessing, title, filename)

## Venn diagram as lots of these gene sets have redundant players

gene_set_list = list(mRNAactivation$gene_symbol,
                     translation_initiation$gene_symbol,
                     translation_elongation$gene_symbol,
                     GCN2AAdeprivation$gene_symbol,
                     translation$gene_symbol,
                     SRPcotranstargeting$gene_symbol,
                     rRNAprocessing$gene_symbol)

category.names = c("mRNA activation",
                   "translation\ninitiation",
                   "translation\nelongation",
                   "GCN2\nAA\ndeprivation",
                   "translation",
                   "SRP\ncotranslational\ntargeting",
                   "rRNA\nprocessing")

ggVennDiagram(
  gene_set_list,
  category.names = category.names,
  label = "count",
  label_geom = "text",
  label_color = "white"
) + scale_fill_gradient(name = "count", trans = "log", na.value = "white")

## plot of the 33 common genes
common33 <- tibble(gene_symbol = Reduce(intersect, gene_set_list)) %>%
  mutate(ensembl_gene = gene_to_ensembl[gene_symbol])

title <- "common 33 genes to all translation related gene sets"
filename <- "02_Output/IC4cyt_translationcommon33.png"
PlotGenesetExpression(cyt_m, common33, title, filename, cluster_cols = T, scale = "none")

## plot of the 11 GCN2 exclusive genes
GetExclusiveGenes <- function(gene_set_list, reference_index){
  ref <- gene_set_list[[reference_index]]
  other <- gene_set_list[-reference_index] %>% flatten_chr() %>% unique()
  exclusive_ref <- setdiff(ref, other)
  return(exclusive_ref)
}

GCN2exclusive <- tibble(gene_symbol = GetExclusiveGenes(gene_set_list, 4)) %>%
  mutate(ensembl_gene = gene_to_ensembl[gene_symbol]) %>% 
  dplyr::filter(gene_symbol != "Gcn1") %>% 
  add_row(gene_symbol = "Gcn1l1", ensembl_gene = gene_to_ensembl["Gcn1l1"]) #this gene changed from Gcn1l1 to Gcn1 in Ensembl anotation and dont find it styarting with gene symbol

title <- "11 genes exclusive to response to GCN2 AA deprivation"
filename <- "02_Output/IC4cyt_GCN2exclusive11.png"
PlotGenesetExpression(cyt_m, GCN2exclusive, title, filename)


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
           show_colnames= FALSE, main = "Stromal transcription of Verginadis et al. 2022 murine CAF markers")

# Samain et al 2021 and Csf1 expression
## Csf1 correlation with ISRact/IC.4.cyt
to_corr <- pivot_longer(cyt_m, cols = !EnsemblID ,names_to = "sample") %>%
  pivot_wider(names_from = EnsemblID, id_cols = sample) %>% 
  dplyr::select(sample, ENSMUSG00000014599) %>%
  inner_join(annotations) %>%
  inner_join(A_cyt)

scatterplot_with_stats <- function(data, varx, vary, label, type, title){
  correlation <- rcorr(data[[varx]], data[[vary]], type = type)
  correlation_txt <- paste0(type, ": R = ", round(correlation$r["x","y"], 2),
                            ", pval = ", round(correlation$P["x","y"], 4))
  ggplot(data) +
    aes_string(x = varx, y = vary, label = label) +
    geom_point() +
    geom_smooth(method=lm) +
    theme_bw() +
    geom_text_repel() +
    labs(title=title, subtitle = correlation_txt)
}

scatterplot_with_stats(to_corr, "ISRact", "ENSMUSG00000014599", "sample", "spearman", "Csf1 mRNA levels as Samain et al. 2021 activation marker")

#scatterplot_with_stats# Houcong et al 2022
### gene weight check
Houcong <- read_tsv("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/01_Input/apCAF_Huocong.tsv", col_names = "apCAF_Huocong")
S_cyt_Huocong_apCAFs <- dplyr::filter(S_cyt, EnsemblID %in% gene_to_ensembl[Houcong$apCAF_Huocong])


ggplot(S_cyt) +
  aes_string("IC.4.cyt") +
  geom_density() +
  geom_vline(xintercept = 0) +
  geom_rug(data = S_cyt_Huocong_apCAFs, aes(IC.4.cyt, color = "Huocong_apCAFs"), outside = T) +
  coord_cartesian(clip = 'off') +
  theme_classic()

## heatmap of gene expression
annot_row <- pivot_longer(Houcong,  cols = tidyselect::everything(), names_to = "signature", values_to = "genes") %>%
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
  dplyr::filter(Genenames %in% flatten_chr(Houcong), !is.na(Genenames)) %>%
  dplyr::select(!EnsemblID) %>%
  relocate(Genenames) %>% column_to_rownames("Genenames") %>% 
  pheatmap(scale = "row", color = colorRampPalette(c("#0C6291", "#FBFEF9", "#A63446"))(100),
           annotation_row = annot_row,
           annotation_col = annot_col,
           annotation_colors = annot_colors,
           cluster_rows = T, show_rownames = TRUE,
           show_colnames= FALSE, main = "Stromal transcription of Verginadis et al. 2022 murine CAF markers")

# Analysis of RNPs 2023
## parameters
omic <- "cyt" # cyt, TEs
interest_IC = "IC.6.TEs"

if (omic == "TEs"){
  expression_m <- TEs_m
} else if (omic == "cyt"){
  expression_m <- cyt_m
}

RBPs <- read_tsv("01_Input/RBPs_mm.csv") %>% # rbp.db mus musculus database Aug 2023
  dplyr::select(-...1) %>%
  dplyr::rename(gene_symbol = `Gene Symbol`, ensembl_gene = `Annotation ID`)

list_of_RBPs <- dplyr::select(RBPs, gene_symbol,ensembl_gene) %>% 
  dplyr::filter(ensembl_gene %in% expression_m$EnsemblID)

## analysis of correlation
RBPs_corrs <-  dplyr::filter(expression_m, EnsemblID %in% list_of_RBPs$ensembl_gene) %>%
  mutate(gene_symbol = ensembl_to_gene[EnsemblID]) %>%
  dplyr::select(gene_symbol, A_TEs$sample) %>% 
  relocate(gene_symbol, all_of(order_by)) %>%
  pivot_longer(-gene_symbol, names_to = "sample") %>%
  pivot_wider(values_from = value, names_from = gene_symbol) %>%
  inner_join(A_TEs, by = "sample") %>%
  column_to_rownames("sample") %>%
  formatted_cors(cor.stat = "spearman") %>%
  filter(measure1 %in% list_of_RBPs$gene_symbol,
         measure2 %in% names(A_TEs)) %>%
  mutate(FDR = p.adjust(p = p, method = "BH"),
         sig_FDR = FDR < 0.05)

write_tsv(RBPs_corrs, paste0("02_Output/RBPs.",omic,"_vs_IC.TEs.tsv")) 
  

## Plots
### n of significant deregulated RBPs per IC
for (i in 1:500){
  random_genes <- sample_n(expression_m, nrow(list_of_RBPs))
  
  random_corr <- random_genes %>%
  dplyr::select(EnsemblID, A_TEs$sample) %>%
  relocate(EnsemblID, all_of(order_by)) %>%
  pivot_longer(-EnsemblID, names_to = "sample") %>%
  pivot_wider(values_from = value, names_from = EnsemblID) %>%
  inner_join(A_TEs, by = "sample") %>%
  column_to_rownames("sample") %>%
  formatted_cors(cor.stat = "spearman") %>%
  filter(measure1 %in% random_genes$EnsemblID,
         measure2 %in% names(A_TEs)) %>%
  mutate(FDR = p.adjust(p = p, method = "BH"),
         sig_FDR = FDR < 0.05)
  
  if (i == 1)
    {
    print(max(random_corr$r))
    random_correlations <- dplyr::mutate(random_corr, iteration = i)
  } else
    {
      print(i)
    random_correlations <- dplyr::mutate(random_corr, iteration = i) %>%
      rbind(random_correlations)
  }
}

random_significant <- dplyr::filter(random_correlations, sig_FDR) %>% 
  group_by(iteration) %>%
  group_by(measure2, .add = TRUE) %>%
  dplyr::summarise(count = n()) %>%
  group_by(measure2) %>%
  dplyr::summarise(mean = mean(count), sd = sd(count))


dplyr::filter(RBPs_corrs, sig_FDR) %>% 
  ggplot(aes(measure2)) +
  geom_bar(stat = "count", fill = "grey69") +
  geom_point(inherit.aes = FALSE, random_significant, mapping = aes(x = measure2, y = mean), color = "red") +
  geom_point(inherit.aes = FALSE, random_significant, mapping = aes(x = measure2, y = mean + sd), shape = 2) +
  geom_point(inherit.aes = FALSE, random_significant, mapping = aes(x = measure2, y = mean - sd), shape = 6) +
  ggtitle(label = paste0("RBPs ", omic, " is significantly correlated to an IC."),
          subtitle = " vs. 500 random lists of genes mean & sd") +
  
  theme_bw() +
  labs(x = "component", y = "count") +
  scale_y_continuous(breaks=seq(0,12,2)) +
  rotate_x_text(angle = 90)

### component specific significantly correlated components
plot_IC6_RBPstotal <- dplyr::filter(RBPs_corrs, sig_FDR,
                                    measure2 == interest_IC) %>%
  arrange(r) %>%
  mutate(measure1 = as_factor(measure1)) %>%
  arrange(abs(r)) %>%
  slice_tail(n=50) %>%
  ggplot(aes(x = r, y = measure1, fill = FDR)) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low = "red", high = "blue", guide=guide_colourbar(reverse = TRUE))+
  theme_bw() +
  xlim(-1,1) +
  labs(title = paste0(interest_IC," top 50 most correlated RBPs ", omic)) +
  xlab("Spearman correlation") +
  ylab("")
  
ggsave(file="02_Output/Figures/plot_IC6TE_RBPstotal.svg", plot=plot_IC6_RBPstotal, width=5, height=5)
