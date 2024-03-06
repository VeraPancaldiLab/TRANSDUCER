#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
library(JADE)
library(Hmisc)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(ggplotify)
library(ggpubr)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/03_IC3Characterization/01_PDX/050324_ISRact_reanalysis/")
source("src/functions.R")

# PARAMETERS
#-------------------------------------------------------------------------------
# Use Remy ICA
use_remy = T 
# Bootstrap #! Not tested for tumour TEs 
run_boot <- FALSE
range.comp <- 2:15 # when ncomp is =< df Warning: In sqrt(puiss[rangeW]) : NaNs produced
boot.iter <- 500
boot.perc <- 0.95

# Analysis
elected_ncomp <- 5 # 4 if looking at distribution, 6 for standar, like in tumour way
component_reorientation = FALSE
reorient <- c(1, 1, -1, 1, 1, -1)
#-------------------------------------------------------------------------------

# DATA LOADING/PROCESSING
#-------------------------------------------------------------------------------
# load TEs
load("../Remy_processed_data/TranslatomeDataForShiny.RData")
TEs <- as_tibble(resi, rownames ="EnsemblID")
sample_info <- read_tsv("~/Documents/02_TRANSDUCER/02_PDX_stroma/00_Data/Processed_data/sample_info.tsv") %>% # Inherited from the stroma analysis
  column_to_rownames("sample") %>% .[colnames(TEs)[-1],]

sample_info$Diabetes <- as_factor(sample_info$Diabetes)

# load annotation with Biomart
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
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

TEs %>% dplyr::filter(EnsemblID %in% non_sex) -> TEs_

## mean centring gene wise
TEs_ %>% column_to_rownames("EnsemblID") %>%
  apply(1, function(x) x - mean(x)) %>% t() %>%
  data.frame() -> TEs__

## Inter quartile range (measure variability)
TEs__ %>% apply(1, IQR) -> iqrs
mostvar <- iqrs[iqrs > median(iqrs)]
TEs__ %>% rownames_to_column("EnsemblID") %>%
  dplyr::filter(EnsemblID %in% names(mostvar)) %>%
  column_to_rownames("EnsemblID") -> TEs_icaready

# ICA Bootstrapping
if (use_remy == F){
  if (run_boot == TRUE){
    ## Baseline
    TEs_icaready %>% jade_range(range.comp, MARGIN = 1) -> base_res_gene
    TEs_icaready %>% jade_range(range.comp, MARGIN = 2) -> base_res_sample
  
    ## Bootstrap
    gene_boot <- jade_choosencom(TEs_icaready, base_res_gene,
                                 MARGIN = 1,
                                 iterations = boot.iter,
                                 seed = 0,
                                 perc = boot.perc
    )
    
    sample_boot <- jade_choosencom(TEs_icaready, base_res_sample,
                                   MARGIN = 2,
                                   iterations = boot.iter,
                                   seed = 0,
                                   perc = boot.perc
    )
    
    boot_plots(s_boot = sample_boot, g_boot = gene_boot, name = "TEs")
  }
  
  # Most robust ICA analysis
  jade_result <- JADE(TEs_icaready, n.comp = elected_ncomp)
  colnames(jade_result[["A"]]) <- paste("IC", 1:elected_ncomp, sep = ".")
  rownames(jade_result[["A"]]) <- names(jade_result$Xmu)
  write_rds(jade_result, "02_Output/tumor_ICA_TEs.RDS") # never save a reoriented ICA analysis


  ## component reorientation
  if (component_reorientation == TRUE){
    stopifnot("vector of reorientation should have the same length as number of components" = length(reorient)==elected_ncomp)
    jade_result[["A"]] <- t(t(jade_result[["A"]]) * reorient)
    jade_result[["S"]] <- t(t(jade_result[["S"]]) * reorient)
  }
  A_mat <- as.data.frame(jade_result[["A"]])
  S_mat <- as.data.frame(jade_result[["S"]])
} else{
  ## Alternatively, use Remy original ICA
  load("../Remy_processed_data/selresICA.RData")
  S_mat = as.data.frame(seltmmICA$S) %>% dplyr::rename_all( ~ paste0("IC.", 1:6)) # manually checked order
  A_mat = as.data.frame(seltmmICA$A) %>% dplyr::rename_all( ~ paste0("IC.", 1:6))
}


#-------------------------------------------------------------------------------

# EXPLORATORY ANALYSIS
#-------------------------------------------------------------------------------

# Sample weight analysis
## Metadata
annotations <- sample_info[-1]
stopifnot(rownames(A_mat) == rownames(annotations))
annotations %>% dplyr::select(!Diabetes) %>%
  names() -> cont_names

## plot
plot_sample_weights(A_mat, annotations, cont_names, "sampleweights_TEs")

## Export for further correlations
stopifnot(rownames(A_mat) == rownames(annotations))
bind_cols(A_mat, annotations[,c("ICA3", "PAMG")]) %>%
  dplyr::rename(ISRact = ICA3) -> complete_annotation

## RNAseq celltype deconvolution
mMCPcounter_res <-  read_tsv("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/02_Output/mMCPcounter_results.tsv") %>% 
  column_to_rownames("cell_types") %>% 
  t() %>% 
  as_tibble(rownames = "samples") %>%
  column_to_rownames("samples") %>% 
  .[rownames(A_mat),]
ImmuCC_res <-  read_tsv("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/02_Output/ImmuCC_results.tsv") %>% column_to_rownames("samples") %>% .[rownames(A_mat),]

Plot_deconv(ImmuCC_res, complete_annotation, "ImmuCC_TEs")
Plot_deconv(mMCPcounter_res, complete_annotation, "mMCPcounter_TEs")

## TF activity
tumour_tf <- read_tsv("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/02_Output/TFact_tumor_PanCan.tsv") %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))
stroma_tf <- read_tsv("~/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/02_Output/TFact_stroma_Gtex.tsv") %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))

Plot_general_TFs(tumour_tf, "TF_tumor_vs_TEs_ICA", 25, complete_annotation)
Plot_general_TFs(stroma_tf, "TF_stroma_vs_TEs_ICA", 25, complete_annotation)
PlotBestCorr(complete_annotation, tumour_tf, 10, "best_TF_tumor_vs_TEs")
PlotBestCorr(complete_annotation, stroma_tf, 10, "best_TF_stroma_vs_TEs")


# Gene weight analysis
## translate to gene names
translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])
reverse_translate = deframe(annot_ensembl75[c("external_gene_id", "ensembl_gene_id")])
## Sanity check

## Sanity Check
Sanity_Check = TEs
Sanity_Check <- dplyr::mutate(Sanity_Check, geneID = translate[EnsemblID]) %>%
  dplyr::filter(!is.na(geneID)) %>% 
  dplyr::select(-EnsemblID) %>%
  distinct(geneID, .keep_all = T) %>%
  column_to_rownames("geneID") %>%
  t()

stopifnot(all(rownames(Sanity_Check) == rownames(annotations)))

Sanity_Check <- cbind(Sanity_Check, complete_annotation) %>%
  rownames_to_column("sampleID")

scatterplot_with_stats(Sanity_Check, varx = "IC.3", vary = "ISRact", type =  "pearson", title = "ISRact vs the object of Remy", label = "sampleID") 
scatterplot_with_stats(Sanity_Check, varx = "IC.3", vary = "ATF4", type =  "pearson", title = "ATF4 vs the object of Remy", label = "sampleID") 

## plot
PlotGeneWeights(S_mat, TEs, 25, translate, complete_annotation, analysis_name = "gene_weights_TEs")

## GSVA of best genes
### gene set preparation
all_genesets <- msigdbr("Homo sapiens")
all_genesets %>% filter(gs_subcat %in% c("CP:REACTOME")) -> use_genesets
msigdbr_list = split(x = use_genesets$ensembl_gene, f = use_genesets$gs_name)

signature_dict <- dplyr::select(all_genesets, c(gs_name, gs_id)) %>%
  distinct(gs_id, .keep_all = T) %>%  deframe()

gsvaRes <- gsva(data.matrix(S_mat), msigdbr_list, min.sz = 15)

### Plot
for (comp in colnames(S_mat)){
  gsvaRes[order(gsvaRes[,comp]),]
  
  gsvaTop <- as_tibble(gsvaRes, rownames = "gene_set") %>% 
    mutate(the_rank = rank(-get(comp), ties.method = "random"),
           gene_set = if_else(str_count(gene_set, "_") < 10, gene_set, signature_dict[gene_set]),
           gene_set = str_remove(gene_set, "REACTOME_"),
           gene_set = fct_reorder(gene_set, the_rank,.desc = T)) %>%
    pivot_longer(cols = -c(gene_set, the_rank), names_to = "component", values_to = "ES") %>% 
    dplyr::filter(the_rank < 25 | the_rank > (nrow(gsvaRes)-25)) %>% 
    mutate(component = if_else(component == comp, comp, "Other")) %>% 
    dplyr::select(!c(the_rank))
  
  ggplot(gsvaTop, aes(x = ES, y = gene_set)) + 
    geom_point(aes(alpha = if_else(component == comp, 0.9, 0.3),
                   color = if_else(ES > 0, "blue", "red"))) +
    theme_bw() +
    labs(title = comp, subtitle = "Best Reactome gene sets") +
    rremove("legend") +
    rremove("ylab")
  ggsave(paste0("02_Output/gsva_TEs", comp, ".pdf"))
}

# # Produce a joint PDF
# concat_pdf <- c("TEs_bootstrapICA_lineplot.pdf",
#                 "TEs_bootstrapICA_boxplot.pdf",
#                 "techdata_corr_TEs.pdf",
#                 "sampleweights_TEs.pdf",
#                 "ImmuCC_TEs.pdf",
#                 "mMCPcounter_TEs.pdf",
#                 "TF_analysis_TF_stroma_vs_TEs_ICA.pdf",
#                 "TF_analysis_TF_tumor_vs_TEs_ICA.pdf")
# 
# for (comp in 1:elected_ncomp){
#   c(concat_pdf, paste("best_TF_stroma_vs_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
#   c(concat_pdf, paste("best_TF_tumor_vs_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
#   c(concat_pdf, paste("gene_weights_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
#   c(concat_pdf, paste("gsva_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
# }
# 
# #system(paste(c("cd 02_Output/ \n pdfunite", concat_pdf, "ICA_TEs.pdf"), collapse = " "))
# system(paste(c("cd 02_Output/ \n convert", concat_pdf, "ICA_TEs.pdf"), collapse = " "))

# RBP expression correlation analysis #! To be done!!
RBPs <- read_tsv("01_Input/RBPs_hs.csv") %>% # rbp.db homo sapiens database March 2024
  dplyr::rename(gene_symbol = `Gene Symbol`, ensembl_gene = `Annotation ID`)

list_of_RBPs <- dplyr::select(RBPs, gene_symbol, ensembl_gene) %>% 
  dplyr::filter(ensembl_gene %in% expression_m$EnsemblID)

## analysis of correlation
# RBPs_corrs <-  dplyr::filter(expression_m, EnsemblID %in% list_of_RBPs$ensembl_gene) %>%
#   mutate(gene_symbol = ensembl_to_gene[EnsemblID]) %>%
#   dplyr::select(gene_symbol, A_TEs$sample) %>% 
#   relocate(gene_symbol, all_of(order_by)) %>%
#   pivot_longer(-gene_symbol, names_to = "sample") %>%
#   pivot_wider(values_from = value, names_from = gene_symbol) %>%
#   inner_join(A_TEs, by = "sample") %>%
#   column_to_rownames("sample") %>%
#   formatted_cors(cor.stat = "spearman") %>%
#   filter(measure1 %in% list_of_RBPs$gene_symbol,
#          measure2 %in% names(A_TEs)) %>%
#   mutate(FDR = p.adjust(p = p, method = "BH"),
#          sig_FDR = FDR < 0.05)
# 
# write_tsv(RBPs_corrs, paste0("02_Output/RBPs.",omic,"_vs_IC.TEs.tsv")) 
# 
# 
# ## Plots
# ### n of significant deregulated RBPs per IC
# for (i in 1:500){
#   random_genes <- sample_n(expression_m, nrow(list_of_RBPs))
#   
#   random_corr <- random_genes %>%
#     dplyr::select(EnsemblID, A_TEs$sample) %>%
#     relocate(EnsemblID, all_of(order_by)) %>%
#     pivot_longer(-EnsemblID, names_to = "sample") %>%
#     pivot_wider(values_from = value, names_from = EnsemblID) %>%
#     inner_join(A_TEs, by = "sample") %>%
#     column_to_rownames("sample") %>%
#     formatted_cors(cor.stat = "spearman") %>%
#     filter(measure1 %in% random_genes$EnsemblID,
#            measure2 %in% names(A_TEs)) %>%
#     mutate(FDR = p.adjust(p = p, method = "BH"),
#            sig_FDR = FDR < 0.05)
#   
#   if (i == 1)
#   {
#     print(max(random_corr$r))
#     random_correlations <- dplyr::mutate(random_corr, iteration = i)
#   } else
#   {
#     print(i)
#     random_correlations <- dplyr::mutate(random_corr, iteration = i) %>%
#       rbind(random_correlations)
#   }
# }
# 
# random_significant <- dplyr::filter(random_correlations, sig_FDR) %>% 
#   group_by(iteration) %>%
#   group_by(measure2, .add = TRUE) %>%
#   dplyr::summarise(count = n()) %>%
#   group_by(measure2) %>%
#   dplyr::summarise(mean = mean(count), sd = sd(count))
# 
# 
# dplyr::filter(RBPs_corrs, sig_FDR) %>% 
#   ggplot(aes(measure2)) +
#   geom_bar(stat = "count") +
#   geom_point(inherit.aes = FALSE, random_significant, mapping = aes(x = measure2, y = mean), color = "red") +
#   geom_point(inherit.aes = FALSE, random_significant, mapping = aes(x = measure2, y = mean + sd), shape = 2) +
#   geom_point(inherit.aes = FALSE, random_significant, mapping = aes(x = measure2, y = mean - sd), shape = 6) +
#   ggtitle(label = paste0("RBPs ", omic, " is significantly correlated to an IC."),
#           subtitle = " vs. 500 random lists of genes mean & sd") +
#   
#   theme_bw() +
#   labs(x = "component", y = "count") +
#   scale_y_continuous(breaks=seq(0,12,2)) +
#   rotate_x_text(angle = 90)
# 
# ### component specific significantly correlated components
# plot_IC6_RBPstotal <- dplyr::filter(RBPs_corrs, sig_FDR,
#                                     measure2 == interest_IC) %>%
#   arrange(r) %>%
#   mutate(measure1 = as_factor(measure1)) %>%
#   arrange(abs(r)) %>%
#   slice_tail(n=50) %>%
#   ggplot(aes(x = r, y = measure1, fill = FDR)) +
#   geom_bar(stat="identity") +
#   scale_fill_gradient(low = "red", high = "blue", guide=guide_colourbar(reverse = TRUE))+
#   theme_bw() +
#   xlim(-1,1) +
#   labs(title = paste0(interest_IC," top 50 most correlated RBPs ", omic)) +
#   xlab("Spearman correlation") +
#   ylab("")
# 
# ggsave(file="02_Output/Figures/plot_IC6TE_RBPstotal.svg", plot=plot_IC6_RBPstotal, width=5, height=5)
# #-------------------------------------------------------------------------------