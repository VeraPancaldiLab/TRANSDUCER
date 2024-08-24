# Package loading
library(tidyverse)
library(mMCPcounter)
# install.packages("devtools")
# library(devtools)
# devtools::install_github("cit-bioinfo/mMCP-counter")
library(Hmisc)
library(corrplot)
library(ggrepel)
setwd("/home/jacobo/Documents/02_TRANSDUCER/02_PDX_stroma/03_Analysis/180122_Various/")

# Data loading
## Expression data
load("../../00_Data/Remy_processed_data/all_RNA/Stroma/Mcpmallrna.RData")
Mcpmallrna # This must be quantile normalized data, probably log2
Mcpmallrna <- Mcpmallrna[,sort(colnames(Mcpmallrna))]
boxplot(Mcpmallrna)

## to gene-IDs
Mcpmallrna_genenames <- Mcpmallrna
ensbl_gene <- read_tsv("01_Input/annotation_ensembl_mm.tsv")
translate <- deframe(ensbl_gene)

Mcpmallrna_genenames <- as_tibble(Mcpmallrna, rownames = "EnsemblID") %>% 
  dplyr::mutate(Gene_ID = translate[EnsemblID]) %>%
  dplyr::select(Gene_ID, matches("PDAC")) %>%
  dplyr::group_by(Gene_ID) %>%
  dplyr::summarise_all(mean) %>%
  dplyr::filter(!is.na(Gene_ID)) %>%
  column_to_rownames("Gene_ID")

## IC3 weights
sample_info <- read_tsv("../../00_Data/Processed_data/sample_info.tsv") %>%
  dplyr::filter(sample %in% colnames(Mcpmallrna)) %>% 
  dplyr::rename(ISRact = "ICA3") %>% 
  column_to_rownames(., var = "sample")

ic3 <- sample_info[, "ISRact", drop =F]

# Deconvolution
## mMCPcounter
?mMCPcounter::mMCPcounter.estimate # runnable with EnsemblIDs
Mcpmallrna_mMCPcounter <- mMCPcounter.estimate(exp = Mcpmallrna,
                                               features = "ENSEMBL.ID",genomeVersion = "GCRm38")
## ImmuCC (CYBERSORT 4 mice)
# download.file(url = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/webserver/SignatureMatrix.rnaseq.csv"
#               destfile = "01_Input/SignatureMatrix.rnaseq.csv")
# unix sed 's/,/\t/g' 01_Input/SignatureMatrix.rnaseq.csv > 01_Input/SignatureMatrix.ImmuCC.tsv
# 
# ### Export to non log .tsv (ensembl IDs as the signature)
# Mcpmallrna_totsv <- as.data.frame(2**Mcpmallrna)
# boxplot(Mcpmallrna_totsv)
# 
# write.table(Mcpmallrna_totsv, "01_Input/Mcpmallrna.tsv",
#             sep = "\t", quote = FALSE, col.names=NA)

Mcpmallrna_ImmuCC <- read_tsv("02_Output/CIBERSORT.Output_Job12.tsv") %>% as.data.frame(.) %>% # For analysis done with "Mcpmallrna.tsv" Job 12. For "HostCyt_foranalysis.tsv" Job 13 (This is used laetr in Deconv)  
  column_to_rownames(., var = "Input Sample") %>% .[,!colnames(.) %in% c("P-value", "Pearson Correlation", "RMSE")]

# Visualization
### mMCPcounter
mcp_plots <- cbind(t(Mcpmallrna_mMCPcounter), as.data.frame(ic3))

#### TCells
mcp_plots.cor <- rcorr(mcp_plots$ISRact, mcp_plots$`T cells`, type = "spearman")
stats <- paste("Spearman: R = ", round(mcp_plots.cor$r["x","y"], 2),
               ", pval = ", round(mcp_plots.cor$P["x","y"], 4), sep = "")

ggplot_mcp_tcells <- ggplot(mcp_plots, aes(x=ISRact, y= `T cells`, label = rownames(mcp_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="MCPcounter", subtitle = stats)
ggsave(file="02_Output/Figures/ggplot_mcp_tcells.svg", plot=ggplot_mcp_tcells, width=5, height=5)

#### CD8 TCells
mcp_plots.cor <- rcorr(mcp_plots$ISRact, mcp_plots$`CD8 T cells`, type = "spearman")
stats <- paste("Spearman: R = ", round(mcp_plots.cor$r["x","y"], 2),
               ", pval = ", round(mcp_plots.cor$P["x","y"], 4), sep = "")

ggplot_mcp_cd8tcells <- ggplot(mcp_plots, aes(x=ISRact, y= `CD8 T cells`, label = rownames(mcp_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="MCPcounter", subtitle = stats)
ggsave(file="02_Output/Figures/ggplot_mcp_cd8tcells.svg", plot=ggplot_mcp_cd8tcells, width=5, height=5)

#### Macrophages
mcp_plots.cor <- rcorr(mcp_plots$ISRact, mcp_plots$`Monocytes / macrophages`, type = "spearman")
stats <- paste("Spearman: R = ", round(mcp_plots.cor$r["x","y"], 2),
               ", pval = ", round(mcp_plots.cor$P["x","y"], 4), sep = "")

ggplot_mcp_monomacro <- ggplot(mcp_plots, aes(x=ISRact, y= `Monocytes / macrophages`, label = rownames(mcp_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="MCPcounter", subtitle = stats)
ggsave(file="02_Output/Figures/ggplot_mcp_monomacro.svg", plot=ggplot_mcp_monomacro, width=5, height=5)



### ImmuCC
immuCC_plots <- cbind(Mcpmallrna_ImmuCC, as.data.frame(ic3))

#### CD4 T Cells
immuCC_plots.cor <-rcorr(immuCC_plots$ISRact, immuCC_plots$`CD4 T Cells`, type = "spearman")
stats <- paste("Spearman: R = ", round(immuCC_plots.cor$r["x","y"], 2),
               ", pval = ", round(immuCC_plots.cor$P["x","y"], 4), sep = "")

ggplot_immucc_cd4tcells <- ggplot(immuCC_plots, aes(x=ISRact, y= `CD4 T Cells`, label = rownames(immuCC_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="immuCC", subtitle = stats)
ggsave(file="02_Output/Figures/ggplot_immucc_cd4tcells.svg", plot=ggplot_immucc_cd4tcells, width=5, height=5)

#### CD8 TCells
immuCC_plots.cor <-rcorr(immuCC_plots$ISRact, immuCC_plots$`CD8 T Cells`, type = "spearman")
stats <- paste("Spearman: R = ", round(immuCC_plots.cor$r["x","y"], 2),
               ", pval = ", round(immuCC_plots.cor$P["x","y"], 4), sep = "")

ggplot_immucc_cd8tcells <- ggplot(immuCC_plots, aes(x=ISRact, y= `CD8 T Cells`, label = rownames(immuCC_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="immuCC", subtitle = stats)
ggsave(file="02_Output/Figures/ggplot_immucc_cd8tcells.svg", plot=ggplot_immucc_cd8tcells, width=5, height=5)


#### Macrophages
immuCC_plots.cor <-rcorr(immuCC_plots$ISRact, immuCC_plots$Macrophages, type = "spearman")
stats <- paste("Spearman: R = ", round(immuCC_plots.cor$r["x","y"], 2),
               ", pval = ", round(immuCC_plots.cor$P["x","y"], 4), sep = "")

ggplot_immucc_macrophages <- ggplot(immuCC_plots, aes(x=ISRact, y= Macrophages, label = rownames(immuCC_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="immuCC", subtitle = stats)
ggsave(file="02_Output/Figures/ggplot_immucc_macrophages.svg", plot=ggplot_immucc_macrophages, width=5, height=5)


### genes
gene_plots <- cbind(t(Mcpmallrna_genenames), as.data.frame(ic3))
gene_plots <- gene_plots[,!duplicated(colnames(gene_plots))]

#### IDO1
gene_plots.cor <-rcorr(gene_plots$ISRact, gene_plots$Ido1, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot_stroma_ido1 <- ggplot(gene_plots, aes(x=ISRact, y=Ido1, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Stroma", subtitle = stats)
ggsave(file="02_Output/Figures/ggplot_stroma_ido1.svg", plot=ggplot_stroma_ido1, width=5, height=5)



#### IL18
gene_plots.cor <-rcorr(gene_plots$ISRact, gene_plots$Il18, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot_stroma_il18 <- ggplot(gene_plots, aes(x=ISRact, y=Il18, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Stroma", subtitle = stats)
ggsave(file="02_Output/Figures/ggplot_stroma_il18.svg", plot=ggplot_stroma_il18, width=5, height=5)



#### IFNG
gene_plots.cor <-rcorr(gene_plots$ISRact, gene_plots$Ifng, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot_stroma_Ifng <- ggplot(gene_plots, aes(x=ISRact, y=Ifng, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Stroma", subtitle = stats)
ggsave(file="02_Output/Figures/ggplot_stroma_Ifng.svg", plot=ggplot_stroma_Ifng, width=5, height=5)


#### PDL1 (CD274)
gene_plots.cor <-rcorr(gene_plots$ISRact, gene_plots$Cd274, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot_stroma_Cd274 <- ggplot(gene_plots, aes(x=ISRact, y=Cd274, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Stroma", subtitle = stats)
ggsave(file="02_Output/Figures/ggplot_stroma_Cd274.svg", plot=ggplot_stroma_Cd274, width=5, height=5)


# TF activity analysis (REDO WITH DATA)
stroma_TFs <- read_tsv("02_Output/TFact_stroma_Gtex.tsv") %>%
  column_to_rownames("TF")

complete_annotation <- dplyr::select(sample_info, PAMG, ISRact, Diabetes)
nTFs = 25
  
correlations <- rcorr(t(stroma_TFs), complete_annotation[,"ISRact",drop=T])
comp_p.adj <- p.adjust(correlations$P[,"y"], "BH")
best_tfs <-  correlations$r[,"y"] %>%
  abs() %>% sort(decreasing = T) %>%
  .[2:(nTFs+1)] %>% names()

annotation_TFs <- tibble(TFs = best_tfs, R = correlations$r[best_tfs,"y"], 
                         p.value = correlations$P[best_tfs,"y"],
                         p.adj = as.numeric(comp_p.adj[best_tfs] < 0.05)) %>% 
  column_to_rownames("TFs")

stroma_TFs %>% .[best_tfs,] %>%
  pheatmap(main="Stroma TFs most correlated with ISRact",
           scale = "row",
           annotation_col = complete_annotation,
           annotation_row = annotation_TFs, 
           filename = "02_Output/Figures/TFact_stroma_ISRact.pdf",
           width = 10,
           height = 8,
           annotation_colors = list(R = c("red", "white", "blue"),
                                    p.value = c("black", "white"),
                                    p.adj = c("white", "black"),
                                    Diabetes = c("white", "black"),
                                    PAMG = c("#FF7F00", "white", "#377DB8"),
                                    ISRact = c("#FFFFCC", "#006837")))

