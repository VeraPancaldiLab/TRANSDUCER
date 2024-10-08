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
library(pdacmolgrad)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/170323_PACAOMICs_ICA/")
source("../100122_ICABoot/functions.R")

# PARAMETERS
#-------------------------------------------------------------------------------
# Bootstrap
run_boot <- FALSE
range.comp <- 2:15 # when ncomp is =< df Warning: In sqrt(puiss[rangeW]) : NaNs produced
boot.iter <- 500
boot.perc <- 0.95

# Analysis
data <- "2017" # "extended" "2017"
elected_ncomp <- 6 # 4 if looking at distribution, 6 for standar, like in tumour way
component_reorientation = TRUE
reorient <- c(1, 1, 1, -1, 1, 1)
#-------------------------------------------------------------------------------

# DATA LOADING/PROCESSING
#-------------------------------------------------------------------------------
# load cyt data
if (data == "2017"){
  PACA_cyt_raw <- read_tsv("../../00_Data/PACAOMICs_data/Nicolle_2017/Murine-Stroma_rawcount_Transcriptome.tsv") %>%
    column_to_rownames("EnsemblID")
  

} else if( data == "extended"){
  PACA_cyt_raw <- read_rds("../../00_Data/PACAOMICs_data/Alexias_53PDX/PDXstromaRaw.rds")
}

## normalize cyt data
PACA_cyt_raw_ <- PACA_cyt_raw[!(apply(PACA_cyt_raw, 1, function(x) {
    any(x == 0)
  })), ] # from anota2seqRemoveZeroSamples()

PACA_cyt <-
  limma::voom(edgeR::calcNormFactors(edgeR::DGEList(PACA_cyt_raw_)))$E %>% ## TMM-log2 from anota2seqNormalize()
  as_tibble(rownames= "EnsemblID")


# load annotation with Biomart for later
## Human (for PAMG)
hs_ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)#listAttributes(ensembl75, page="feature_page")

hs_annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id'), mart = hs_ensembl75)

hs_translate <- deframe(hs_annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

## mice
ensembl75 <- useEnsembl(biomart = "genes",
                        dataset = "mmusculus_gene_ensembl",
                        version = 75)#listAttributes(ensembl75, page="feature_page")

annot_ensembl75 <- getBM(attributes = c('ensembl_gene_id',
                                        'external_gene_id',
                                        'entrezgene',
                                        'mgi_id',
                                        'chromosome_name'), mart = ensembl75)

translate <- deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

## Create sample_info_PACAOMICS
sample_info_ <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>%
  dplyr::select(- RNAconc_cyt, -RNAconc_pol, -cyt_tumor_counts, -pol_tumor_counts, -cyt_host_counts, -pol_host_counts, -PAMG)

sample_info_$Diabetes <- as_factor(sample_info_$Diabetes)

### add True PAMG
if (data == "2017"){
  tumor_counts <- read_tsv("../../00_Data/PACAOMICs_data/Nicolle_2017/Human-Tumor_rawcount_Transcriptome.tsv") %>%
    mutate(GeneNames = hs_translate[EnsemblID]) %>% 
    dplyr::select(-EnsemblID) %>% 
    group_by(GeneNames) %>% 
    summarise_all(mean) %>%
    column_to_rownames("GeneNames") %>%
    DGEList() %>%
    calcNormFactors(method= "upperquartile") %>%
    cpm(log=TRUE)

}else if (data == "extended"){
  load("../../00_Data/PACAOMICs_data/Alexias_53PDX/PDX_HUMAN_RAW.RData")

  tumor_counts <- as_tibble(x, rownames = "EnsemblID") %>% 
    mutate(GeneNames = hs_translate[EnsemblID]) %>% 
    dplyr::select(-EnsemblID) %>% 
    group_by(GeneNames) %>% 
    summarise_all(mean) %>%
    column_to_rownames("GeneNames") %>%
    DGEList() %>%
    calcNormFactors(method= "upperquartile") %>%
    cpm(log=TRUE)
}

type_pamg <- projectMolGrad(newexp = tumor_counts,
                            geneSymbols = rownames(tumor_counts)) %>%
  as_tibble(rownames = "sample") %>%
  dplyr::select(sample, PDX) %>%
  dplyr::rename(PAMG = PDX)

sample_info <- right_join(sample_info_, type_pamg, by = "sample") %>% column_to_rownames("sample")
  

# filtering
## XY
annot_ensembl75 %>% dplyr::filter(!chromosome_name %in% c("X", "Y")) %>%
  pull(ensembl_gene_id) -> non_sex

PACA_cyt %>% dplyr::filter(EnsemblID %in% non_sex) -> PACA_cyt_

## mean centring gene wise
PACA_cyt_ %>% column_to_rownames("EnsemblID") %>%
  apply(1, function(x) x - mean(x)) %>% t() %>%
  data.frame() -> PACA_cyt__

## Inter quartile range (measure variability)
PACA_cyt__ %>% apply(1, IQR) -> iqrs
mostvar <- iqrs[iqrs > median(iqrs)]
PACA_cyt__ %>% rownames_to_column("EnsemblID") %>%
  dplyr::filter(EnsemblID %in% names(mostvar)) %>%
  column_to_rownames("EnsemblID") -> PACA_cyt_icaready

# ICA Bootstrapping
if (run_boot == TRUE){
  ## Baseline
  PACA_cyt_icaready %>% jade_range(range.comp, MARGIN = 1) -> base_res_gene
  PACA_cyt_icaready %>% jade_range(range.comp, MARGIN = 2) -> base_res_sample
  
  ## Bootstrap
  gene_boot <- jade_choosencom(PACA_cyt_icaready, base_res_gene,
                               MARGIN = 1,
                               iterations = boot.iter,
                               seed = 0,
                               perc = boot.perc
  )
  
  sample_boot <- jade_choosencom(PACA_cyt_icaready, base_res_sample,
                                 MARGIN = 2,
                                 iterations = boot.iter,
                                 seed = 0,
                                 perc = boot.perc
  )
  
  boot_plots(s_boot = sample_boot, g_boot = gene_boot, name = "PACA_cyt")
} 

# Most robust ICA analysis
jade_result <- JADE(PACA_cyt_icaready, n.comp = elected_ncomp)
colnames(jade_result[["A"]]) <- paste("IC", 1:elected_ncomp, sep = ".")
rownames(jade_result[["A"]]) <- names(jade_result$Xmu)
write_rds(jade_result, "02_Output/ICA_PACA_cyt.RDS") # never save a reoriented ICA analysis

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
annotations <- sample_info[-1] %>% .[rownames(A_mat),]
stopifnot(rownames(A_mat) == rownames(annotations))
annotations %>% dplyr::select(!Diabetes) %>%
  names() -> cont_names

## plot
plot_sample_weights(A_mat = A_mat, annotations = annotations, cont_names = cont_names, analysis_name = "sampleweights_PACA_cyt")


# FIXED UNTIL HERE!!!


## Export for further correlations
stopifnot(rownames(A_mat) == rownames(annotations))
complete_annotation <- bind_cols(A_mat, annotations[,c("ICA3", "PAMG")]) %>%
  dplyr::rename(ISRact = ICA3)

## RNAseq celltype deconvolution
mMCPcounter_res <-  read_tsv("../180122_Various/02_Output/mMCPcounter_results.tsv") %>% column_to_rownames("cell_types") %>% t() %>% as_tibble(rownames = "samples") %>% column_to_rownames("samples") %>% .[rownames(A_mat),]
ImmuCC_res <-  read_tsv("../180122_Various/02_Output/ImmuCC_results.tsv") %>% column_to_rownames("samples") %>% .[rownames(A_mat),]

Plot_deconv(ImmuCC_res, complete_annotation, "ImmuCC_PACA_cyt")
Plot_deconv(mMCPcounter_res, complete_annotation, "mMCPcounter_PACA_cyt")

## TF activity
tumour_tf <- read_tsv("../180122_Various/02_Output/TFact_tumor_PanCan.tsv") %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))
stroma_tf <- read_tsv("../180122_Various/02_Output/TFact_stroma_Gtex.tsv") %>% column_to_rownames("TF") %>% dplyr::select(rownames(complete_annotation))

Plot_general_TFs(tumour_tf, "TF_tumor_vs_PACA_cyt_ICA", 25, complete_annotation)
Plot_general_TFs(stroma_tf, "TF_stroma_vs_PACA_cyt_ICA", 25, complete_annotation)
PlotBestCorr(complete_annotation, tumour_tf, 10, "best_TF_tumor_vs_PACA_cyt")
PlotBestCorr(complete_annotation, stroma_tf, 10, "best_TF_stroma_vs_PACA_cyt")

# Gene weight analysis
## plot
PlotGeneWeights(S_mat, PACA_cyt, 25, translate, complete_annotation, analysis_name = "gene_weights_PACA_cyt")

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
  ggsave(paste0("02_Output/gsva_PACA_cyt", comp, ".pdf"))
  
  # gsvaTop %>% pheatmap(filename = paste("02_Output/gsva_PACA_cyt", comp, ".pdf", sep=""), height = 10 , width = 15,
  #                      cluster_rows = F, main = paste(comp, "\n Best gene sets")) # to have a pheatmap
  # gsvaTop  %>% rownames() %>% msigdb_descriptions[.,] %>%
  #   cbind(gsvaTop[comp], .) %>% rownames_to_column("msigdb") %>% write_tsv("02_Output/gsvaRes_IC11nc13.tsv") #to write
}

# Network analysis
MID <- read_csv("01_Input/MID2022.csv")

trans_mgi_name <- deframe(annot_ensembl75[c("mgi_id", "external_gene_id")])
trans_ensembl_name <- deframe(annot_ensembl75[c("ensembl_gene_id", "external_gene_id")])

nodes <- as_tibble(S_mat, rownames = "ensembl_id") %>%
  mutate(id = trans_ensembl_name[ensembl_id]) %>%
  relocate(ensembl_id, .after = "id")

edges <- mutate(MID,
                gene1 = trans_mgi_name[gene1],
                gene2 = trans_mgi_name[gene2]) %>%
  dplyr::filter(gene1 %in% nodes$id, 
                gene2 %in% nodes$id) %>%
  dplyr::rename(source = gene1,
                target = gene2,
                interaction = type) %>%
  dplyr::select(source, target, interaction)

PlotNetwork(nodes, edges, S_mat,valid_comp = 1:6, main_name = "MID_ICAPACA_cyt")

#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: MCPCounter deconvolution # CHECK IF THIS WORKS
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
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAPACA_cyt ~ MCPCounter",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/deconvolution_corrplot.svg", plot=deconvolution_corrplot, width=8, height=6)
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
  labs(x = NULL, y = NULL, fill = "Spearman's\nAbsolute\nCorrelation", title="Correlations ICAPACA_cyt ~ MCPCounter",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
  scale_y_discrete(expand=c(0,0), limits = paste("IC", rev(1:elected_ncomp), sep = ".")) +
  ggpubr::rotate_x_text(angle = 90)

ggsave(file="02_Output/Figures/stroma_tf_corrplot.svg", plot=stroma_tf_corrplot, width=10, height=6)

#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: all the correlation plots of above together
#-------------------------------------------------------------------------------
sampleweight_corrplots <- ggarrange(plotlist = list(deconvolution_corrplot, stroma_tf_corrplot), ncol = 2, common.legend = T, align = "h", widths = c(0.4, 0.5), legend = "right")
ggsave(file="02_Output/Figures/sampleweight_corrplots.svg", plot=sampleweight_corrplots, width=18, height=6)

#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: IC.4 Gene weights bottom and top genes
#-------------------------------------------------------------------------------
# Get best genes
n_genes = 10
S_mat %>% arrange(desc(IC.4)) -> S_sort
S_sort %>% tail(10) -> IC4_mneg
S_sort %>% head(10) -> IC4_mpos

# Density plot
S_mat %>% ggplot() +
  aes_string(y = "IC.4") +
  geom_density(alpha=.5, fill="lightgrey") +
  scale_x_reverse() +
  scale_y_continuous(expand = c(0,0))+
  geom_hline(yintercept = 0, colour = "black") +
  annotate("rect",xmin = -Inf, xmax = Inf,   ymin = min(IC4_mneg["IC.4"]), ymax = max(IC4_mneg["IC.4"]),   fill = "blue", alpha = 0.5) +
  annotate("rect",xmin = -Inf, xmax = Inf,   ymin =  min(IC4_mpos["IC.4"]), ymax = max(IC4_mpos["IC.4"]),   fill = "red", alpha = 0.5) +
  theme_classic() -> densplot

# Heatmap
## Annotation
annot_row_ <- tibble(name = rownames(IC4_mpos), class = "most possitive")
tibble(name = rownames(IC4_mneg), class = "most negative") %>% bind_rows(annot_row_) %>% column_to_rownames("name") -> annot_row__
annot_row <- annot_row__
rownames(annot_row) <- rownames(annot_row__) %>% translate[.] %>% make.names(unique = TRUE)
  
annot_col <- complete_annotation[c("PAMG", "ISRact", "IC.4")]

annot_colors <- list(class = c(`most possitive` = "red", `most negative` = "blue"),
                     PAMG = c("#FF7F00", "white", "#377DB8"),
                     ISRact = c("#FFFFCC", "#006837"),
                     IC.4 = c("#FFFFCC", "#5b0066")) 

## Get gene names
ensembl_toplot_ <- PACA_cyt %>% dplyr::filter(EnsemblID %in% c(rownames(IC4_mpos), rownames(IC4_mneg))) %>%
  arrange(match(EnsemblID, c(rownames(IC4_mpos), rownames(IC4_mneg))))

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
density_heatmap_IC4 <- ggarrange(densplot + rremove("xylab"), heatmap, heights = c(1.5, 10), widths = c(0.2,1),
                                 ncol = 2, nrow = 1)

ggsave(file="02_Output/Figures/density_heatmap_IC4.svg", plot=density_heatmap_IC4, width=10, height=6)

#-------------------------------------------------------------------------------

# FIGURE SPECIFIC PLOTS: IC.4 GSVA bottom and top genes
#-------------------------------------------------------------------------------
gsvaRes[order(gsvaRes[,"IC.4"]),]

gsvaTop <- as_tibble(gsvaRes, rownames = "gene_set") %>% 
  mutate(the_rank = rank(-IC.4, ties.method = "random"),
         #gene_set = if_else(str_count(gene_set, "_") < 10, gene_set, signature_dict[gene_set]),
         gene_set = str_remove(gene_set, "REACTOME_"),
         gene_set = fct_reorder(gene_set, the_rank,.desc = T)) %>%
  pivot_longer(cols = -c(gene_set, the_rank), names_to = "component", values_to = "ES") %>% 
  dplyr::filter(the_rank < 15 | the_rank > (nrow(gsvaRes)-15)) %>% 
  mutate(component = if_else(component == "IC.4", "IC.4", "Other")) %>% 
  dplyr::select(!c(the_rank))

gsva_IC4 <- ggplot(gsvaTop, aes(x = ES, y = gene_set)) + 
  geom_point(aes(alpha = if_else(component == "IC.4", 0.9, 0.3),
                 color = if_else(ES > 0, "blue", "red"))) +
  theme_bw() +
  labs(title = "IC.4", subtitle = "Best Reactome gene sets") +
  rremove("legend") +
  rremove("ylab")

ggsave(file="02_Output/Figures/gsva_IC4.svg", plot=gsva_IC4, width=10, height=6)
