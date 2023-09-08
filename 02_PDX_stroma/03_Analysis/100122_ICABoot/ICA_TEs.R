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

# PARAMETERS
#-------------------------------------------------------------------------------
# Bootstrap
run_boot <- FALSE
range.comp <- 2:15 # when ncomp is =< df Warning: In sqrt(puiss[rangeW]) : NaNs produced
boot.iter <- 500
boot.perc <- 0.95

# Analysis
elected_ncomp <- 6 # 4 if looking at distribution, 6 for standar, like in tumour way
component_reorientation = TRUE
reorient <- c(1, 1, -1, 1, 1, -1)
#-------------------------------------------------------------------------------

# DATA LOADING/PROCESSING
#-------------------------------------------------------------------------------
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
                          'mgi_id',
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
write_rds(jade_result, "02_Output/ICA_TEs.RDS") # never save a reoriented ICA analysis

## component reorientation
if (component_reorientation == TRUE){
  stopifnot("vector of reorientation should have the same length as number of components" = length(reorient)==elected_ncomp)
  jade_result[["A"]] <- t(t(jade_result[["A"]]) * reorient)
  jade_result[["S"]] <- t(t(jade_result[["S"]]) * reorient)
}
A_mat <- as.data.frame(jade_result[["A"]])
S_mat <- as.data.frame(jade_result[["S"]])
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

## QC and alignment info
multiqc <- read_tsv("../../00_Data/Processed_data/multiQC_summary.tsv") %>%
  dplyr::filter(CITID %in% names(TEs)[-1])


seq_info <- read_tsv("../../00_Data/Processed_data/sequencing_info.tsv") %>%
  dplyr::filter(CITID %in% names(TEs)[-1])

## Cyt read info
multiqc_cyt <- dplyr::filter(multiqc, Fraction == "Cytosolic") %>%
  dplyr::select(-c(fastq_name, Fraction))

multiqc_cyt_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(multiqc_cyt) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(multiqc_cyt), measure2 %in% names(A_mat)) %>% #not square corr
  mutate(r = abs(r)) %>% #abscorr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICATEs ~ multiqc",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggpubr::rotate_x_text(angle = 90)


seq_info_cyt <- dplyr::filter(seq_info, Fraction == "Cytosolic") %>%
  mutate(input2maped_length_diff = Average_input_read_length - Average_mapped_length) %>%
  dplyr::select(-c(fastq_id, Fraction, Average_input_read_length, Average_mapped_length))

star_cyt_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(seq_info_cyt) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(seq_info_cyt), measure2 %in% names(A_mat)) %>% #not square corr
  mutate(r = abs(r)) %>% #abscorr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICATEs ~ STAR metadata",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggpubr::rotate_x_text(angle = 90)

techi_plots_cyt <- ggarrange(multiqc_cyt_corr, star_cyt_corr, align = "h", common.legend = T, legend = "right")

## Polysome (F8-9) read info
multiqc_pol <- dplyr::filter(multiqc, Fraction == "Polysome") %>%
  dplyr::select(-c(fastq_name, Fraction))

multiqc_pol_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(multiqc_pol) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(multiqc_pol), measure2 %in% names(A_mat)) %>% #not square corr
  mutate(r = abs(r)) %>% #abscorr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICATEs ~ multiqc",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggpubr::rotate_x_text(angle = 90)

seq_info_pol <- dplyr::filter(seq_info, Fraction == "Polysome") %>%
  mutate(input2maped_length_diff = Average_input_read_length - Average_mapped_length) %>%
  dplyr::select(-c(fastq_id, Fraction, Average_input_read_length, Average_mapped_length))

star_pol_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(seq_info_pol) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(seq_info_pol), measure2 %in% names(A_mat)) %>% #not square corr
  mutate(r = abs(r)) %>% #abscorr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICATEs ~ STAR metadata",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(0,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  ggpubr::rotate_x_text(angle = 90)

techi_plots_pol <- ggarrange(multiqc_pol_corr, star_pol_corr, align = "h", common.legend = T, legend = "right")

## Export
techi_plots <- ggarrange(techi_plots_cyt, techi_plots_pol,
                         labels = c("Cytosolic", "Polysome"),
                         font.label = list(face = "bold", color = "black", size = 20),
                         vjust=-1,  nrow = 2, align = "v")

annotate_figure(techi_plots,
                top = text_grob("", size = 40))

ggsave("02_Output/techdata_corr_TEs.pdf", dpi = "print", width = 210, height = 300, units = "mm", scale = 1.25)

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

concat_pdf <- c("TEs_bootstrapICA_lineplot.pdf",
                "TEs_bootstrapICA_boxplot.pdf",
                "techdata_corr_TEs.pdf",
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

PlotNetwork(nodes, edges, S_mat,valid_comp = 1:6, main_name = "MID_ICATE")

#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: Clinical and technical data
#-------------------------------------------------------------------------------
columns_multiqc = c("CITID", "Dups_R1", "GC_R1")
multiqc_sp <- dplyr::inner_join(dplyr::select(multiqc_cyt, all_of(columns_multiqc)),
                                dplyr::select(multiqc_pol, all_of(columns_multiqc)),
                                by="CITID", suffix = c("_cyt", "_pol"))

clinical_technical_corrplot <- complete_annotation %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(multiqc_sp) %>%
  dplyr::select(CITID, matches("IC."), ISRact, PAMG, matches("Dups_R1"), matches("GC_R1")) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% c("ISRact", "PAMG", "Dups_R1_cyt", "GC_R1_cyt", "Dups_R1_pol", "GC_R1_pol"),
         measure2 %in% names(A_mat)) %>% #not square corr
  ggplot(aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAcyt ~ multiqc",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/clinical_technical_corrplot_TEs.svg", plot=clinical_technical_corrplot, width=4, height=6)
sed -i "s/ textLength='[^']*'//" file.svg
#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: MCPCounter deconvolution
#-------------------------------------------------------------------------------
deconvolution_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(as_tibble(mMCPcounter_res, rownames = "CITID")) %>%
  dplyr::select(CITID, names(A_mat), names(mMCPcounter_res)) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(mMCPcounter_res),
         measure2 %in% names(A_mat))

clust <- dplyr::select(deconvolution_corr, measure1, measure2, r) %>% 
  pivot_wider(names_from=measure2, values_from = r) %>%
  column_to_rownames("measure1") %>%
  dist() %>%
  hclust()

deconvolution_corrplot <- ggplot(deconvolution_corr, aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAcyt ~ MCPCounter",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/deconvolution_corrplot_TEs.svg", plot=deconvolution_corrplot, width=8, height=6)
#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: Most correlated stromal TF activities
#-------------------------------------------------------------------------------
stroma_mostvar_TF <- Get_mostvar(stroma_tf, 20) %>% rownames_to_column() %>% pivot_longer(-rowname) %>% 
  pivot_wider(names_from=rowname, values_from=value)

stroma_tf_corr <- A_mat %>%
  as_tibble(rownames = "CITID") %>%
  inner_join(as_tibble(t(stroma_tf), rownames = "CITID")) %>%
  dplyr::select(CITID, names(A_mat), names(stroma_mostvar_TF)[-1]) %>%
  column_to_rownames("CITID") %>%
  formatted_cors("spearman") %>%
  filter(measure1 %in% names(stroma_mostvar_TF)[-1],
         measure2 %in% names(A_mat))

clust <- dplyr::select(stroma_tf_corr, measure1, measure2, r) %>% 
  pivot_wider(names_from=measure2, values_from = r) %>%
  column_to_rownames("measure1") %>%
  dist() %>%
  hclust()

stroma_tf_corrplot <- ggplot(stroma_tf_corr, aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAcyt ~ MCPCounter",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/stroma_tf_corrplot_TEs.svg", plot=stroma_tf_corrplot, width=10, height=6)

#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: all the correlation plots of above together
#-------------------------------------------------------------------------------
sampleweight_corrplots <- ggarrange(plotlist = list(clinical_technical_corrplot, deconvolution_corrplot, stroma_tf_corrplot), ncol = 3, common.legend = T, align = "h", widths = c(0.165, 0.4, 0.5), legend = "right")
ggsave(file="02_Output/Figures/sampleweight_corrplots_TEs.svg", plot=sampleweight_corrplots, width=20, height=6)

#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: IC.6 Gene weights bottom and top genes
#-------------------------------------------------------------------------------
# Get best genes
n_genes = 10
S_mat %>% arrange(desc(IC.6)) -> S_sort
S_sort %>% tail(10) -> IC6_mneg
S_sort %>% head(10) -> IC6_mpos

# Density plot
S_mat %>% ggplot() +
  aes_string(y = "IC.6") +
  geom_density(alpha=.5, fill="lightgrey") +
  scale_x_reverse() +
  scale_y_continuous(expand = c(0,0))+
  geom_hline(yintercept = 0, colour = "black") +
  annotate("rect",xmin = -Inf, xmax = Inf,   ymin = min(IC6_mneg["IC.6"]), ymax = max(IC6_mneg["IC.6"]),   fill = "blue", alpha = 0.5) +
  annotate("rect",xmin = -Inf, xmax = Inf,   ymin =  min(IC6_mpos["IC.6"]), ymax = max(IC6_mpos["IC.6"]),   fill = "red", alpha = 0.5) +
  theme_classic() -> densplot

# Heatmap
## Annotation
annot_row_ <- tibble(name = rownames(IC6_mpos), class = "most possitive")
tibble(name = rownames(IC6_mneg), class = "most negative") %>% bind_rows(annot_row_) %>% column_to_rownames("name") -> annot_row__
annot_row <- annot_row__
rownames(annot_row) <- rownames(annot_row__) %>% translate[.] %>% make.names(unique = TRUE)

annot_col <- complete_annotation[c("PAMG", "ISRact", "IC.6")]

annot_colors <- list(class = c(`most possitive` = "red", `most negative` = "blue"),
                     PAMG = c("#FF7F00", "white", "#377DB8"),
                     ISRact = c("#FFFFCC", "#006837"),
                     IC.6 = c("#FFFFCC", "#5b0066")) 

## Get gene names
ensembl_toplot_ <- TEs %>% dplyr::filter(EnsemblID %in% c(rownames(IC6_mpos), rownames(IC6_mneg))) %>%
  arrange(match(EnsemblID, c(rownames(IC6_mpos), rownames(IC6_mneg))))

genes_toplot <- ensembl_toplot_
genes_toplot$Genenames <- ensembl_toplot_$EnsemblID %>% translate[.] %>%
  make.names(unique = TRUE)

## Plot
heatmap <- genes_toplot %>% dplyr::select(!EnsemblID) %>%
  relocate(Genenames) %>% column_to_rownames("Genenames") %>% 
  pheatmap(scale = "row", color = colorRampPalette(c("#0C6291", "#FBFEF9", "#A63446"))(100),
           annotation_row = annot_row, annotation_col = annot_col, annotation_colors = annot_colors,
           cluster_rows = F, show_colnames = TRUE) %>%
  as.grob()

# Export joint figure
density_heatmap_IC6 <- ggarrange(densplot + rremove("xylab"), heatmap, heights = c(1.5, 10), widths = c(0.2,1),
                                 ncol = 2, nrow = 1)

ggsave(file="02_Output/Figures/density_heatmap_IC6_TEs.svg", plot=density_heatmap_IC6, width=10, height=6)
