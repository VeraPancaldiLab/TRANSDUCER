#!/usr/bin/env Rscript
library(tidyverse)
library(biomaRt)
library(JADE)
library(corrplot)
library(Hmisc)
library(pheatmap)
library(msigdbr)
library(GSVA)
library(ggplotify)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/")
source("functions.R")
run_boot <- TRUE
range.comp <- 2:15 # when ncomp is =< df Warning: In sqrt(puiss[rangeW]) : NaNs produced
boot.iter <- 500

# load TEs
TEs <- read_tsv("../081221_TranslationEfficacy/02_Output/TEs.tsv")
sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>%
  column_to_rownames("sample") %>% .[colnames(TEs)[-1],]

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
if (run_boot == TRUE){
  ## Baseline
  TEs_icaready %>% jade_range(range.comp, MARGIN = 1) -> base_res_gene
  TEs_icaready %>% jade_range(range.comp, MARGIN = 2) -> base_res_sample

  ## Bootstrap
  gene_boot <- jade_choosencom(TEs_icaready, base_res_gene,
                               MARGIN = 1,
                               iterations = boot.iter,
                               seed = 0
  )
  
  sample_boot <- jade_choosencom(TEs_icaready, base_res_sample,
                                 MARGIN = 2,
                                 iterations = boot.iter,
                                 seed = 0
  )
  
  boot_plots(s_boot = sample_boot, g_boot = gene_boot, name = "TEs")
}

# Most robust ICA analysis
elected_ncomp <- 9 

jade_result <- JADE(TEs_icaready, n.comp = elected_ncomp)
colnames(jade_result[["A"]]) <- paste("IC", 1:elected_ncomp, sep = ".")
rownames(jade_result[["A"]]) <- names(jade_result$Xmu)

# Sample weight analysis
## Metadata
A_mat <- as.data.frame(jade_result[["A"]])
S_mat <- as.data.frame(jade_result[["S"]])
annotations <- sample_info[-1]
stopifnot(rownames(A_mat) == rownames(annotations))
annotations %>% dplyr::select(!Diabetes) %>%
  names() -> cont_names

## plot
plot_sample_weights(A_mat, annotations, cont_names, "sampleweights_TEs")

## Export for further correlations
stopifnot(rownames(A_mat) == rownames(annotations))
bind_cols(A_mat, annotations[,c("ICA3", "PAMG")]) %>%
  rename(SerDep = ICA3) -> complete_annotation

## RNAseq celltype deconvolution
mMCPcounter_res <-  read_tsv("../180122_Various/02_Output/mMCPcounter_results.tsv") %>% column_to_rownames("cell_types") %>% t() %>% as_tibble(rownames = "samples") %>% column_to_rownames("samples") %>% .[rownames(A_mat),]
ImmuCC_res <-  read_tsv("../180122_Various/02_Output/ImmuCC_results.tsv") %>% column_to_rownames("samples") %>% .[rownames(A_mat),]

Plot_deconv(ImmuCC_res, complete_annotation, "ImmuCC_TEs")
Plot_deconv(mMCPcounter_res, complete_annotation, "mMCPcounter_TEs")

## TF activity
tumour_tf <- read_tsv("../180122_Various/02_Output/TFact_tumor_PanCan.tsv") %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))
stroma_tf <- read_tsv("../180122_Various/02_Output/TFact_stroma_Gtex.tsv") %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))

Plot_general_TFs(tumour_tf, "TF_tumor_vs_TEs_ICA", 25, complete_annotation)
Plot_general_TFs(stroma_tf, "TF_stroma_vs_TEs_ICA", 25, complete_annotation)
PlotBestCorr(complete_annotation, tumour_tf, 10, "best_TF_tumor_vs_TEs")
PlotBestCorr(complete_annotation, stroma_tf, 10, "best_TF_stroma_vs_TEs")

# Gene weight analysis
## translate to gene names
translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

## plot
PlotGeneWeights(S_mat, TEs, 25, translate, complete_annotation, analysis_name = "gene_weights_TEs")

## GSVA of best genes
### gene set preparation
all_genesets <- msigdbr("Mus musculus")
all_genesets %>% filter(gs_cat %in% c("H", "C2", "C7")) -> use_genesets
msigdbr_list = split(x = use_genesets$ensembl_gene, f = use_genesets$gs_name)
msigdb_descriptions <- use_genesets[c("gs_name", "gs_description")] %>%
  unique() %>% column_to_rownames("gs_name")

gsvaRes <- gsva(data.matrix(S_mat), msigdbr_list, min.sz = 15)

### Plot
for (comp in colnames(S_mat)){
  gsvaRes[order(gsvaRes[,comp]),]
  
  gsvaRes %>% as.data.frame() %>%  mutate(the_rank = rank(-get(comp), ties.method = "random")) %>% 
    dplyr::filter(the_rank < 25 | the_rank > (nrow(gsvaRes)-25)) %>% arrange(the_rank) %>%
    dplyr::select(!c(the_rank)) -> gsvaTop
  
  
  gsvaTop %>% pheatmap(filename = paste("02_Output/gsva_TEs", comp, ".pdf", sep=""), height = 10 , width = 15,
                       cluster_rows = F, main = paste(comp, "\n Best gene sets")) # to have a pheatmap
  # gsvaTop  %>% rownames() %>% msigdb_descriptions[.,] %>% 
  #   cbind(gsvaTop[comp], .) %>% rownames_to_column("msigdb") %>% write_tsv("02_Output/gsvaRes_IC11nc13.tsv") #to write
}


concat_pdf <- c("TEs_bootstrapICA_lineplot.pdf",
                "TEs_bootstrapICA_boxplot.pdf",
                "sampleweights_TEs.pdf",
                "ImmuCC_TEs.pdf",
                "mMCPcounter_TEs.pdf",
                "TF_analysis_TF_stroma_vs_TEs_ICA.pdf",
                "TF_analysis_TF_tumor_vs_TEs_ICA.pdf")

for (comp in 1:elected_ncomp){
  c(concat_pdf, paste("best_TF_stroma_vs_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
  c(concat_pdf, paste("best_TF_tumor_vs_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
  c(concat_pdf, paste("gene_weights_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
  c(concat_pdf, paste("gsva_TEsIC.", comp, ".pdf", sep = "")) -> concat_pdf
}

#system(paste(c("cd 02_Output/ \n pdfunite", concat_pdf, "ICA_TEs.pdf"), collapse = " "))
system(paste(c("cd 02_Output/ \n convert", concat_pdf, "ICA_TEs.pdf"), collapse = " "))
