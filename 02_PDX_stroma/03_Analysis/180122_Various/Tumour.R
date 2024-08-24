# Package loading
library(tidyverse)
# install.packages("devtools")
# library(devtools)
# devtools::install_github("cit-bioinfo/mMCP-counter")
library(Hmisc)
library(corrplot)
library(ggrepel)
library(pheatmap)
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/")

# Data loading
## Expression data
load("../../00_Data/Remy_processed_data/all_RNA/Tumeur/Hcpmallrna.RData")

Hcpmallrna # This must be quantile normalized data, probably log2
Hcpmallrna <- Hcpmallrna[,sort(colnames(Hcpmallrna))]
boxplot(Hcpmallrna)

## to gene-IDs
Hcpmallrna_genenames <- Hcpmallrna
ensbl_gene <- read_tsv("01_Input/annotation_ensembl_hs.tsv")
translate <- deframe(ensbl_gene)

Hcpmallrna %>% rownames(.) %>% translate[.] -> rownames(Hcpmallrna_genenames)
Hcpmallrna_genenames <- Hcpmallrna_genenames[!is.na(rownames(Hcpmallrna_genenames)),]

## IC3 weights
sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>%
  dplyr::filter(sample %in% colnames(Hcpmallrna)) %>% 
  dplyr::rename(ISRact = "ICA3") %>% 
  column_to_rownames(., var = "sample")

ic3 <- sample_info[, "ISRact", drop =F]


### genes
all(rownames(t(Hcpmallrna_genenames)) == rownames(ic3))
gene_plots <- cbind(t(Hcpmallrna_genenames), as.data.frame(ic3))
gene_plots <- gene_plots[,!duplicated(colnames(gene_plots))] # needed to plot

#### IL18
gene_plots.cor <-rcorr(gene_plots$ISRact, gene_plots$IL18, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot_tumour_il18 <- ggplot(gene_plots, aes(x=ISRact, y=IL18, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor", subtitle = stats)

ggsave(file="02_Output/Figures/ggplot_tumour_il18.svg", plot=ggplot_tumour_il18, width=5, height=5)

#### IDO1
gene_plots.cor <- rcorr(gene_plots$ISRact, gene_plots$IDO1, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot_tumour_IDO1 <- ggplot(gene_plots, aes(x=ISRact, y=IDO1, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor", subtitle = stats)

ggsave(file="02_Output/Figures/ggplot_tumour_IDO1.svg", plot=ggplot_tumour_IDO1, width=5, height=5)

#### CD274
gene_plots.cor <- rcorr(gene_plots$ISRact, gene_plots$CD274, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot_tumour_CD274 <- ggplot(gene_plots, aes(x=ISRact, y=CD274, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor", subtitle = stats)

ggsave(file="02_Output/Figures/ggplot_tumour_CD274.svg", plot=ggplot_tumour_CD274, width=5, height=5)

# TF activity analysis (REDO WITH DATA)
tumour_TFs <- read_tsv("02_Output/TFact_tumor_PanCan.tsv") %>%
  column_to_rownames("TF")

complete_annotation <- dplyr::select(sample_info, PAMG, ISRact, Diabetes)
nTFs = 25

correlations <- rcorr(t(tumour_TFs), complete_annotation[,"ISRact",drop=T])
comp_p.adj <- p.adjust(correlations$P[,"y"], "BH")
best_tfs <-  correlations$r[,"y"] %>%
  abs() %>% sort(decreasing = T) %>%
  .[2:(nTFs+1)] %>% names()

annotation_TFs <- tibble(TFs = best_tfs, R = correlations$r[best_tfs,"y"], 
                         p.value = correlations$P[best_tfs,"y"],
                         p.adj = as.numeric(comp_p.adj[best_tfs] < 0.05)) %>% 
  column_to_rownames("TFs")

tumour_TFs %>% .[best_tfs,] %>%
  pheatmap(main="Tumour TFs most correlated with ISRact",
           scale = "row",
           annotation_col = complete_annotation,
           annotation_row = annotation_TFs, 
           filename = "02_Output/Figures/TFact_tumour_ISRact.pdf",
           width = 10,
           height = 8,
           annotation_colors = list(R = c("red", "white", "blue"),
                                    p.value = c("black", "white"),
                                    p.adj = c("white", "black"),
                                    Diabetes = c("white", "black"),
                                    PAMG = c("#FF7F00", "white", "#377DB8"),
                                    ISRact = c("#FFFFCC", "#006837")))

#-------------------------------------------------------------------------------
#5vs5 most extreme samples TF activities
subset_annotation <- dplyr::arrange(complete_annotation, desc(ISRact)) %>%
  rownames_to_column("sampleid") %>%
  rowid_to_column() %>%
  dplyr::mutate(group = if_else(rowid<6, "high_ISRact", if_else(rowid>(27-5), "low_ISRact", "intermediate_ISRact"))) %>%
  dplyr::filter(group != "intermediate_ISRact") %>%
  dplyr::select(-rowid) %>%
  column_to_rownames("sampleid")

subset_TFs <- tumour_TFs[rownames(subset_annotation)]
nTFs = 25
correlations <- rcorr(t(subset_TFs), subset_annotation[,"ISRact",drop=T])
comp_p.adj <- p.adjust(correlations$P[,"y"], "BH")
best_tfs <-  correlations$r[,"y"] %>%
  abs() %>% sort(decreasing = T) %>%
  .[2:(nTFs+1)] %>% names()

subset_TFs %>% .[best_tfs,]  %>%
  pheatmap(main="Tumour TFs most correlated with ISRact",
           scale = "row",
           annotation_col = subset_annotation,
           annotation_row = annotation_TFs, 
           filename = "02_Output/Figures/TFact_tumour_5vs5_ISRact.pdf",
           width = 7,
           height = 5,
           annotation_colors = list(R = c("red", "white", "blue"),
                                    p.value = c("black", "white"),
                                    p.adj = c("white", "black"),
                                    Diabetes = c("white", "black"),
                                    PAMG = c("#FF7F00", "white", "#377DB8"),
                                    ISRact = c("#FFFFCC", "#006837"),
                                    group = c(low_ISRact = "seagreen", high_ISRact= "tomato3")))

#-------------------------------------------------------------------------------