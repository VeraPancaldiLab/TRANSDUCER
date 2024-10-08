# Package loading
library(tidyverse)
library(mMCPcounter)
# install.packages("devtools")
# library(devtools)
# devtools::install_github("cit-bioinfo/mMCP-counter")
library(Hmisc)
library(corrplot)
library(ggrepel)
library(GSVA)
setwd("/home/jacobo/Documents/02_TRANSDUCER/03_IC3Characterization/01_PDX/270521_Stroma_IFNG/")

# Data loading
## Expression data
load("../Remy_processed_data/all_RNA/Stroma/Mcpmallrna.RData")
Mcpmallrna # This must be quantile normalized data, probably log2
Mcpmallrna <- Mcpmallrna[,sort(colnames(Mcpmallrna))]
boxplot(Mcpmallrna)

## to gene-IDs
Mcpmallrna_genenames <- Mcpmallrna
ensbl_gene <- read_tsv("01_Input/GRCm39_ensemblvsgenename.txt")
translate <- deframe(ensbl_gene)

Mcpmallrna %>% rownames(.) %>% translate[.] -> rownames(Mcpmallrna_genenames)
Mcpmallrna_genenames <- Mcpmallrna_genenames[!is.na(rownames(Mcpmallrna_genenames)),]

## IC3 weights
ic3 <- read_tsv("../Remy_processed_data/samplesIC3_custom.csv") %>%
  as.data.frame(.) %>%
  column_to_rownames(., var = "CITID") %>%
  .[colnames(Mcpmallrna),"ICA3SampleWeight", drop = FALSE]


ic3

# Deconvolution
## mMCPcounter devtools::install_github("cit-bioinfo/mMCP-counter")
?mMCPcounter::mMCPcounter.estimate # runnable with EnsemblIDs
Mcpmallrna_mMCPcounter <- mMCPcounter.estimate(exp = Mcpmallrna,
                                               features = "ENSEMBL.ID")
## ImmuCC (CYBERSORT 4 mice)
# download.file(url = "https://raw.githubusercontent.com/wuaipinglab/ImmuCC/master/webserver/SignatureMatrix.rnaseq.csv"
#               destfile = "01_Input/SignatureMatrix.rnaseq.csv")
#unix sed 's/,/\t/g' 01_Input/SignatureMatrix.rnaseq.csv > 01_Input/SignatureMatrix.ImmuCC.tsv

### Export to non log .tsv (ensembl IDs as the signature)
Mcpmallrna_totsv <- as.data.frame(2**Mcpmallrna)
boxplot(Mcpmallrna_totsv)

write.table(Mcpmallrna_totsv, "01_Input/Mcpmallrna.tsv",
            sep = "\t", quote = FALSE, col.names=NA)

Mcpmallrna_ImmuCC <- read_tsv("02_Output/CIBERSORT.Output_Job12.tsv") %>% as.data.frame(.) %>%
  column_to_rownames(., var = "Input Sample") %>% .[,!colnames(.) %in% c("P-value", "Pearson Correlation", "RMSE")]


## Visualization
### correlation ic3 - mMCPcunter
mmcp_ic3 <- rcorr(t(Mcpmallrna_mMCPcounter), data.matrix(ic3), type = "spearman")

corrplot(mmcp_ic3$r, order="hclust", type = "lower", p.mat = mmcp_ic3$P,
         sig.level = 0.05, method="color", insig = "label_sig",
         title = "Stroma Composition (mMCPcounter) vs IC3", mar=c(0,0,1,0))# http://stackoverflow.com/a/14754408/54964

### correlation ic3 - ImmuCC
immucc_ic3 <- rcorr(data.matrix(Mcpmallrna_ImmuCC), data.matrix(ic3), type = "spearman")

corrplot(immucc_ic3$r, order="hclust", type = "lower", p.mat = immucc_ic3$P,
         sig.level = 0.05, method="color", insig = "label_sig",
         title = "Stroma Composition vs IC3 (ImmuCC~CIBERSORT)", mar=c(0,0,1,0))# http://stackoverflow.com/a/14754408/54964

### correlation mMCPcunter - ImmuCC
mmcp_immuCC <- rcorr(t(Mcpmallrna_mMCPcounter), data.matrix(Mcpmallrna_ImmuCC), type = "spearman")
mmcp_immuCC$r <- mmcp_immuCC$r[rownames(Mcpmallrna_mMCPcounter), colnames(Mcpmallrna_ImmuCC)]
mmcp_immuCC$P <- mmcp_immuCC$P[rownames(Mcpmallrna_mMCPcounter), colnames(Mcpmallrna_ImmuCC)]

corrplot(mmcp_immuCC$r, order="original", type = "full", p.mat = mmcp_immuCC$P,
         sig.level = 0.05, method="color", insig = "label_sig",
         title = "Stroma Composition mMCPc. vs ImmuCC", mar=c(0,0,1,0))# http://stackoverflow.com/a/14754408/54964

### scaterplots
#### mMCPcounter
mcp_plots <- cbind(t(Mcpmallrna_mMCPcounter), as.data.frame(ic3))

##### TCells
mcp_plots.cor <-rcorr(mcp_plots$ICA3SampleWeight, mcp_plots$`T cells`, type = "spearman")
stats <- paste("Spearman: R = ", round(mcp_plots.cor$r["x","y"], 2),
               ", pval = ", round(mcp_plots.cor$P["x","y"], 4), sep = "")

ggplot(mcp_plots, aes(x=ICA3SampleWeight, y= `T cells`, label = rownames(mcp_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="MCPcounter", subtitle = stats)

##### CD8 TCells
mcp_plots.cor <-rcorr(mcp_plots$ICA3SampleWeight, mcp_plots$`CD8 T cells`, type = "spearman")
stats <- paste("Spearman: R = ", round(mcp_plots.cor$r["x","y"], 2),
               ", pval = ", round(mcp_plots.cor$P["x","y"], 4), sep = "")

ggplot(mcp_plots, aes(x=ICA3SampleWeight, y= `CD8 T cells`, label = rownames(mcp_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="MCPcounter", subtitle = stats)



##### Macrophages
mcp_plots.cor <-rcorr(mcp_plots$ICA3SampleWeight, mcp_plots$`Monocytes / macrophages`, type = "spearman")
stats <- paste("Spearman: R = ", round(mcp_plots.cor$r["x","y"], 2),
               ", pval = ", round(mcp_plots.cor$P["x","y"], 4), sep = "")

ggplot(mcp_plots, aes(x=ICA3SampleWeight, y= `Monocytes / macrophages`, label = rownames(mcp_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="MCPcounter", subtitle = stats)


#### ImmuCC
immuCC_plots <- cbind(Mcpmallrna_ImmuCC, as.data.frame(ic3))

##### CD4 T Cells
immuCC_plots.cor <-rcorr(immuCC_plots$ICA3SampleWeight, immuCC_plots$`CD4 T Cells`, type = "spearman")
stats <- paste("Spearman: R = ", round(immuCC_plots.cor$r["x","y"], 2),
               ", pval = ", round(immuCC_plots.cor$P["x","y"], 4), sep = "")

ggplot(immuCC_plots, aes(x=ICA3SampleWeight, y= `CD4 T Cells`, label = rownames(immuCC_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="immuCC", subtitle = stats)

##### CD8 TCells
immuCC_plots.cor <-rcorr(immuCC_plots$ICA3SampleWeight, immuCC_plots$`CD8 T Cells`, type = "spearman")
stats <- paste("Spearman: R = ", round(immuCC_plots.cor$r["x","y"], 2),
               ", pval = ", round(immuCC_plots.cor$P["x","y"], 4), sep = "")

ggplot(immuCC_plots, aes(x=ICA3SampleWeight, y= `CD8 T Cells`, label = rownames(immuCC_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="immuCC", subtitle = stats)


##### Macrophages
immuCC_plots.cor <-rcorr(immuCC_plots$ICA3SampleWeight, immuCC_plots$`Macrophages`, type = "spearman")
stats <- paste("Spearman: R = ", round(immuCC_plots.cor$r["x","y"], 2),
               ", pval = ", round(immuCC_plots.cor$P["x","y"], 4), sep = "")

ggplot(immuCC_plots, aes(x=ICA3SampleWeight, y= `Macrophages`, label = rownames(immuCC_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="immuCC", subtitle = stats)

# Check presence of neurons
## Axon markers by GSEA of Farias et al 2020 bioinfo analysis 
##of common 300 axon specific genes (https://rnajournal.cshlp.org/content/26/5/595.full)
axon_gs <- read_tsv("01_Input/Farias2020_300axongenes.tsv") %>%
  rename_with(~c("EnsemblID", "Description", "GeneName"), everything())

axon_enrich <- gsva(gsvaParam(exprData = Mcpmallrna, geneSets = as.list(axon_gs["EnsemblID"]), minSize = 15))["EnsemblID",]

## Visualization
### scatterplot
axonb_plots <- as_tibble(ic3, rownames = "samples") %>% mutate(Farias2020_coreaxon = axon_enrich[samples])

axonb_plots.cor <-rcorr(axonb_plots$ICA3SampleWeight, axonb_plots$Farias2020_coreaxon, type = "spearman")
stats <- paste("Spearman: R = ", round(axonb_plots.cor$r["x","y"], 2),
               ", pval = ", round(axonb_plots.cor$P["x","y"], 4), sep = "")

ggplot(axonb_plots, aes(x=ICA3SampleWeight, y=Farias2020_coreaxon, label = axonb_plots$samples)) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Stroma", subtitle = stats)

### heatmap
Mcpmallrna[axon_gs$EnsemblID,] %>% 
  pheatmap(scale = "row", annotation_col = column_to_rownames(axonb_plots, "samples"),
           show_rownames = F)
  
# Individual gene plots
### genes
gene_plots <- cbind(t(Mcpmallrna_genenames), as.data.frame(ic3))
gene_plots <- gene_plots[,!duplicated(colnames(gene_plots))]

#### PHGDH
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$Phgdh, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=Phgdh, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Stroma", subtitle = stats)

#### CBS
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$Cbs, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=Cbs, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Stroma", subtitle = stats)

#### PSAT1
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$Psat1, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=Psat1, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Stroma", subtitle = stats)

#### IDO1
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$Ido1, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=Ido1, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Stroma", subtitle = stats)

