# Package loading

library(tidyverse)
library(mMCPcounter)
# install.packages("devtools")
# library(devtools)
# devtools::install_github("cit-bioinfo/mMCP-counter")
library(Hmisc)
library(corrplot)
library(ggrepel)
setwd("/home/jacobo/Documents/03_Sauyeun_paper/01_PDX/270521_Stroma_IFNG/")

# Data loading
## Expression data
load("../Remy_processed_data/all_RNA/Tumeur/Hcpmallrna.RData")

Hcpmallrna # This must be quantile normalized data, probably log2
Hcpmallrna <- Hcpmallrna[,sort(colnames(Hcpmallrna))]
boxplot(Hcpmallrna)

## to gene-IDs
Hcpmallrna_genenames <- Hcpmallrna
ensbl_gene <- read_tsv("01_Input/GRCh38.p13_ensemblvsgenename.txt")
translate <- deframe(ensbl_gene)

Hcpmallrna %>% rownames(.) %>% translate[.] -> rownames(Hcpmallrna_genenames)
Hcpmallrna_genenames <- Hcpmallrna_genenames[!is.na(rownames(Hcpmallrna_genenames)),]

## IC3 weights
ic3 <- read_tsv("../Remy_processed_data/samplesIC3_custom.csv") %>%
  as.data.frame(.) %>%
  column_to_rownames(., var = "CITID") %>%
  .[colnames(Hcpmallrna),"ICA3SampleWeight", drop = FALSE]


### genes
all(rownames(t(Hcpmallrna_genenames)) == rownames(ic3))
gene_plots <- cbind(t(Hcpmallrna_genenames), as.data.frame(ic3))
gene_plots <- gene_plots[,!duplicated(colnames(gene_plots))] # needed to plot

#### PHGDH
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$PHGDH, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=PHGDH, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor", subtitle = stats)

#### CBS
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$CBS, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=CBS, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor", subtitle = stats)

#### PSAT1
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$PSAT1, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=PSAT1, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor", subtitle = stats)


#### IDO1
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$IDO1, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=IDO1, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor", subtitle = stats)


#### IL18
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$IL18, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=IL18, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor", subtitle = stats)


#### IFNG
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$IFNG, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=IFNG, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor", subtitle = stats)


#### PDL1 (AK CD274)
gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, gene_plots$CD274, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=ICA3SampleWeight, y=CD274, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor", subtitle = stats)


# Find IC3 indicator
## PHGDH & CBS
ggplot(gene_plots, aes(x=CBS, y=PHGDH, color= ICA3SampleWeight, label = rownames(gene_plots))) + 
  geom_point() + 
#  scale_x_continuous(trans='log10') +
#  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "PHGDH/CBS as IC3 indicator")

## Only CBS
ggplot(gene_plots, aes(x=CBS, y=ICA3SampleWeight, color= ICA3SampleWeight, label = rownames(gene_plots))) + 
  geom_point() + 
#  scale_x_continuous(trans='log10') +
#  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "CBS as IC3 indicator")


## Only PHGDH
ggplot(gene_plots, aes(x=PHGDH, y=ICA3SampleWeight, color= ICA3SampleWeight, label = rownames(gene_plots))) + 
  geom_point() + 
  #  scale_x_continuous(trans='log10') +
  #  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "PHGDH as IC3 indicator")


## Only IL18
ggplot(gene_plots, aes(x=IL18, y=ICA3SampleWeight, color= ICA3SampleWeight, label = rownames(gene_plots))) + 
  geom_point() + 
  #  scale_x_continuous(trans='log10') +
  #  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "IL18 as IC3 indicator")

## IL18 and PHGDH
ggplot(gene_plots, aes(x=IL18, y=PHGDH, color= ICA3SampleWeight, label = rownames(gene_plots))) + 
  geom_point() + 
  #  scale_x_continuous(trans='log10') +
  #  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "IL18/PHGDH as IC3 indicator")

## combined IL18PHGDH
### Addition
gene_plots$comIL18PHGDH1 <- scale(gene_plots$IL18) - scale(gene_plots$PHGDH)

gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, c(gene_plots$comIL18PHGDH1), type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")
print(stats)
ggplot(gene_plots, aes(x=comIL18PHGDH1, y=ICA3SampleWeight, color= ICA3SampleWeight, label = rownames(gene_plots))) + 
  geom_point() + 
  #  scale_x_continuous(trans='log10') +
  #  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "combined IL18 PHGDH (IL18 - PHGDH) as IC3 indicator")


### Mean
gene_plots$comIL18PHGDH2 <- rowMeans(cbind(gene_plots$IL18,-gene_plots$PHGDH))

gene_plots.cor <-rcorr(gene_plots$ICA3SampleWeight, c(gene_plots$comIL18PHGDH2), type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")
print(stats)
ggplot(gene_plots, aes(x=comIL18PHGDH2, y=ICA3SampleWeight, color= ICA3SampleWeight, label = rownames(gene_plots))) + 
  geom_point() + 
  #  scale_x_continuous(trans='log10') +
  #  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "combined IL18 PHGDH mean(IL18, -PHGDH) as IC3 indicator")

