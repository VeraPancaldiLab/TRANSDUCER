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
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/100122_ICABoot/")
source("functions.R")
run_boot <- FALSE
range.comp <- 2:15 # when ncomp is =< df Warning: In sqrt(puiss[rangeW]) : NaNs produced
boot.iter <- 500
boot.perc <- 0.95

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
                                        'mgi_id',
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

# ICA Bootstrapping
if (run_boot == TRUE){
  ## Baseline
  cyt_icaready %>% jade_range(range.comp, MARGIN = 1) -> base_res_gene
  cyt_icaready %>% jade_range(range.comp, MARGIN = 2) -> base_res_sample
  
  ## Bootstrap
  gene_boot <- jade_choosencom(cyt_icaready, base_res_gene,
                               MARGIN = 1,
                               iterations = boot.iter,
                               seed = 0,
                               perc = boot.perc
  )
  
  sample_boot <- jade_choosencom(cyt_icaready, base_res_sample,
                                 MARGIN = 2,
                                 iterations = boot.iter,
                                 seed = 0,
                                 perc = boot.perc
  )
  
  boot_plots(s_boot = sample_boot, g_boot = gene_boot, name = "Cyt")
} 

# Most robust ICA analysis
elected_ncomp <- 6 # 4 if looking at distribution, 6 for standar, like in tumour way

jade_result <- JADE(cyt_icaready, n.comp = elected_ncomp)
colnames(jade_result[["A"]]) <- paste("IC", 1:elected_ncomp, sep = ".")
rownames(jade_result[["A"]]) <- names(jade_result$Xmu)
write_rds(jade_result, "02_Output/ICA_cyt.RDS")

# Sample weight analysis
## Metadata
A_mat <- as.data.frame(jade_result[["A"]])
S_mat <- as.data.frame(jade_result[["S"]])
annotations <- sample_info[-1]
stopifnot(rownames(A_mat) == rownames(annotations))
annotations %>% dplyr::select(!Diabetes) %>%
  names() -> cont_names

## plot
plot_sample_weights(A_mat = A_mat, annotations = annotations, cont_names = cont_names, analysis_name = "sampleweights_cyt")

## Export for further correlations
stopifnot(rownames(A_mat) == rownames(annotations))
bind_cols(A_mat, annotations[,c("ICA3", "PAMG")]) %>%
  dplyr::rename(ISRact = ICA3) -> complete_annotation

## QC and alignment info
multiqc <- read_tsv("../../00_Data/Processed_data/multiQC_summary.tsv") %>%
  dplyr::filter(CITID %in% names(cyt)[-1], Fraction == "Cytosolic") %>%
  dplyr::select(-c(fastq_name, Fraction))

multiqc_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(multiqc) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(multiqc), measure2 %in% names(A_mat)) %>% #not square corr
  mutate(r = abs(r)) %>% #abscorr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAcyt ~ multiqc",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggpubr::rotate_x_text(angle = 90)


seq_info <- read_tsv("../../00_Data/Processed_data/sequencing_info.tsv") %>%
  dplyr::filter(CITID %in% names(cyt)[-1], Fraction == "Cytosolic") %>%
  mutate(input2maped_length_diff = Average_input_read_length - Average_mapped_length) %>%
  dplyr::select(-c(fastq_id, Fraction, Average_input_read_length, Average_mapped_length))

star_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(seq_info) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(seq_info), measure2 %in% names(A_mat)) %>% #not square corr
  mutate(r = abs(r)) %>% #abscorr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAcyt ~ STAR metadata",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggpubr::rotate_x_text(angle = 90)

techi_plots <- ggarrange(multiqc_corr, star_corr, align = "h", common.legend = T, legend = "right")
annotate_figure(techi_plots, top = text_grob("Cytosolic info",face = "bold", color = "black", size = 20))
ggsave("02_Output/techdata_corr_cyt.pdf", dpi = "print", width = 210, height = 150, units = "mm", scale = 1.25)

## RNAseq celltype deconvolution
mMCPcounter_res <-  read_tsv("../180122_Various/02_Output/mMCPcounter_results.tsv") %>% column_to_rownames("cell_types") %>% t() %>% as_tibble(rownames = "samples") %>% column_to_rownames("samples") %>% .[rownames(A_mat),]
ImmuCC_res <-  read_tsv("../180122_Various/02_Output/ImmuCC_results.tsv") %>% column_to_rownames("samples") %>% .[rownames(A_mat),]

Plot_deconv(ImmuCC_res, complete_annotation, "ImmuCC_cyt")
Plot_deconv(mMCPcounter_res, complete_annotation, "mMCPcounter_cyt")

## TF activity
tumour_tf <- read_tsv("../180122_Various/02_Output/TFact_tumor_PanCan.tsv") %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))
stroma_tf <- read_tsv("../180122_Various/02_Output/TFact_stroma_Gtex.tsv") %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))

Plot_general_TFs(tumour_tf, "TF_tumor_vs_cyt_ICA", 25, complete_annotation)
Plot_general_TFs(stroma_tf, "TF_stroma_vs_cyt_ICA", 25, complete_annotation)
PlotBestCorr(complete_annotation, tumour_tf, 10, "best_TF_tumor_vs_cyt")
PlotBestCorr(complete_annotation, stroma_tf, 10, "best_TF_stroma_vs_cyt")

# Gene weight analysis
## translate to gene names
translate = deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

## plot
PlotGeneWeights(S_mat, cyt, 25, translate, complete_annotation, analysis_name = "gene_weights_cyt")

## GSVA of best genes
### gene set preparation
all_genesets <- msigdbr("Mus musculus")
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
           gene_set = fct_reorder(gene_set, the_rank)) %>%
    pivot_longer(cols = -c(gene_set, the_rank), names_to = "component", values_to = "ES") %>% 
    dplyr::filter(the_rank < 25 | the_rank > (nrow(gsvaRes)-25)) %>% 
    mutate(component = if_else(component == comp, comp, "Other")) %>% 
    dplyr::select(!c(the_rank))

  ggplot(gsvaTop, aes(x = ES, y = gene_set)) + 
    geom_point(aes(alpha = if_else(component == comp, 0.9, 0.3),
                   color = if_else(ES > 0, "red", "blue"))) +
    theme_bw() +
    labs(title = comp, subtitle = "Best Reactome gene sets") +
    rremove("legend") +
    rremove("ylab")
  ggsave(paste0("02_Output/gsva_cyt", comp, ".pdf"))
  
  # gsvaTop %>% pheatmap(filename = paste("02_Output/gsva_cyt", comp, ".pdf", sep=""), height = 10 , width = 15,
  #                      cluster_rows = F, main = paste(comp, "\n Best gene sets")) # to have a pheatmap
  # gsvaTop  %>% rownames() %>% msigdb_descriptions[.,] %>%
  #   cbind(gsvaTop[comp], .) %>% rownames_to_column("msigdb") %>% write_tsv("02_Output/gsvaRes_IC11nc13.tsv") #to write
}


concat_pdf <- c("Cyt_bootstrapICA_lineplot.pdf",
                "Cyt_bootstrapICA_boxplot.pdf",
                "techdata_corr_cyt.pdf",
                "sampleweights_cyt.pdf",
                "ImmuCC_cyt.pdf",
                "mMCPcounter_cyt.pdf",
                "TF_analysis_TF_stroma_vs_cyt_ICA.pdf",
                "TF_analysis_TF_tumor_vs_cyt_ICA.pdf")

for (comp in 1:elected_ncomp){
  c(concat_pdf, paste("best_TF_stroma_vs_cytIC.", comp, ".pdf", sep = "")) -> concat_pdf
  c(concat_pdf, paste("best_TF_tumor_vs_cytIC.", comp, ".pdf", sep = "")) -> concat_pdf
  c(concat_pdf, paste("gene_weights_cytIC.", comp, ".pdf", sep = "")) -> concat_pdf
  c(concat_pdf, paste("gsva_cytIC.", comp, ".pdf", sep = "")) -> concat_pdf
}

#system(paste(c("cd 02_Output/ \n pdfunite", concat_pdf, "ICA_cyt.pdf"), collapse = " "))
system(paste(c("cd 02_Output/ \n convert", concat_pdf, "ICA_cyt.pdf"), collapse = " "))


# Network analysis
MID <- read_csv("01_Input/MID2022.csv")

trans_mgi_ensembl <- deframe(annot_ensembl75[c("mgi_id", "ensembl_gene_id")])
trans_ensembl_name <- deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

nodes <- as_tibble(S_mat, rownames = "id") %>%
  mutate(gene_name = trans_ensembl_name[id]) %>%
  relocate(gene_name, .after = "id")

edges <- mutate(MID,
                gene1 = trans_mgi_ensembl[gene1],
                gene2 = trans_mgi_ensembl[gene2]) %>%
  dplyr::filter(gene1 %in% nodes$id, 
                gene2 %in% nodes$id) %>%
  dplyr::rename(source = gene1,
                target = gene2,
                interaction = type) %>%
  dplyr::select(source, target, interaction)

PlotNetwork(nodes, edges, S_mat,valid_comp = 1:6, main_name = "MID_ICAcyt")
