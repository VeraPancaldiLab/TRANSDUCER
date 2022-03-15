# Package loading
library(tidyverse)
library(mMCPcounter)
# install.packages("devtools")
# library(devtools)
# devtools::install_github("cit-bioinfo/mMCP-counter")
library(Hmisc)
library(corrplot)
library(ggrepel)
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/03_IC3Characterization/01_PDX/270521_Stroma_IFNG/")

# Data loading
## Expression data (cyt)
load("../Remy_processed_data/all_RNA/Tumeur/Hcpmallrna.RData")

Hcpmallrna # This must be quantile normalized data, probably log2
Hcpmallrna <- Hcpmallrna[,sort(colnames(Hcpmallrna))]
boxplot(Hcpmallrna)

## Polysome (from shiny app)
load("../Remy_processed_data/TranslatomeDataForShiny.RData") # shiny app data (Tumour residuals included)
ribo %>% as_tibble(rownames = "EnsemblID") %>%
  column_to_rownames("EnsemblID") %>% as.matrix() -> pol

## TEs (from shiny app)
resi %>% as_tibble(rownames = "EnsemblID") %>%
  column_to_rownames("EnsemblID") %>% as.matrix() -> TEs

## to gene-IDs
ensbl_gene <- read_tsv("01_Input/GRCh38.p13_ensemblvsgenename.txt")
translate <- deframe(ensbl_gene)

Hcpmallrna_genenames <- Hcpmallrna
Hcpmallrna %>% rownames(.) %>% translate[.] -> rownames(Hcpmallrna_genenames)
Hcpmallrna_genenames <- Hcpmallrna_genenames[!is.na(rownames(Hcpmallrna_genenames)),]

pol_genenames <- pol
pol %>% rownames(.) %>% translate[.] -> rownames(pol_genenames)
pol_genenames <- pol_genenames[!is.na(rownames(pol_genenames)),]

TEs_genenames <- TEs
TEs %>% rownames(.) %>% translate[.] -> rownames(TEs_genenames)
TEs_genenames <- TEs_genenames[!is.na(rownames(TEs_genenames)),]

## IC3 weights
read_tsv("../Remy_processed_data/samplesIC3_custom.csv") %>%
  dplyr::select(CITID,ICA3SampleWeight) %>% dplyr::rename(SerDep = ICA3SampleWeight) %>%
  dplyr::arrange(match(CITID, colnames(Hcpmallrna))) %>%
  column_to_rownames("CITID") -> ic3

### preparation to plot
all(rownames(t(Hcpmallrna_genenames)) == rownames(ic3))
gene_plots <- cbind(t(Hcpmallrna_genenames), as.data.frame(ic3))
gene_plots <- gene_plots[,!duplicated(colnames(gene_plots))] # needed to plot

all(rownames(t(pol_genenames)) == rownames(ic3))
pol_plots <- cbind(t(pol_genenames), as.data.frame(ic3))
pol_plots <- pol_plots[,!duplicated(colnames(pol_plots))]

all(rownames(t(TEs_genenames)) == rownames(ic3))
TEs_plots <- cbind(t(TEs_genenames), as.data.frame(ic3))
TEs_plots <- TEs_plots[,!duplicated(colnames(TEs_plots))]

### Cyt expression
#### PHGDH
gene_plots.cor <-rcorr(gene_plots$SerDep, gene_plots$PHGDH, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=SerDep, y=PHGDH, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor cyt", subtitle = stats)

#### CBS
gene_plots.cor <-rcorr(gene_plots$SerDep, gene_plots$CBS, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=SerDep, y=CBS, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor cyt", subtitle = stats)

#### PSAT1
gene_plots.cor <-rcorr(gene_plots$SerDep, gene_plots$PSAT1, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=SerDep, y=PSAT1, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor cyt", subtitle = stats)


#### IDO1
gene_plots.cor <-rcorr(gene_plots$SerDep, gene_plots$IDO1, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=SerDep, y=IDO1, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor cyt", subtitle = stats)


#### IL18
gene_plots.cor <-rcorr(gene_plots$SerDep, gene_plots$IL18, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=SerDep, y=IL18, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor cyt", subtitle = stats)


#### IFNG
gene_plots.cor <-rcorr(gene_plots$SerDep, gene_plots$IFNG, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=SerDep, y=IFNG, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor cyt", subtitle = stats)


#### PDL1 (AK CD274)
gene_plots.cor <-rcorr(gene_plots$SerDep, gene_plots$CD274, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=SerDep, y=CD274, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor cyt", subtitle = stats)

### DUSP6 
gene_plots.cor <-rcorr(gene_plots$SerDep, gene_plots$DUSP6, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=SerDep, y=DUSP6, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor cyt", subtitle = stats)

#### CD73 AKA NT5E (immunomodulator)
gene_plots.cor <-rcorr(gene_plots$SerDep, gene_plots$NT5E, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=SerDep, y=NT5E, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor cyt", subtitle = stats)

#### JOSD2 (degrade PHGDH)
gene_plots.cor <-rcorr(gene_plots$SerDep, gene_plots$JOSD2, type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")

ggplot(gene_plots, aes(x=SerDep, y=JOSD2, label = rownames(gene_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor cyt", subtitle = stats)


## Polysome fraction
### ATF4 (sanity check)
pol_plots.cor <-rcorr(pol_plots$SerDep, pol_plots$ATF4, type = "spearman")
stats <- paste("Spearman: R = ", round(pol_plots.cor$r["x","y"], 2),
               ", pval = ", round(pol_plots.cor$P["x","y"], 4), sep = "")

ggplot(pol_plots, aes(x=SerDep, y=ATF4, label = rownames(pol_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor pol", subtitle = stats)

### DUSP6 
pol_plots.cor <-rcorr(pol_plots$SerDep, pol_plots$DUSP6, type = "spearman")
stats <- paste("Spearman: R = ", round(pol_plots.cor$r["x","y"], 2),
               ", pval = ", round(pol_plots.cor$P["x","y"], 4), sep = "")

ggplot(pol_plots, aes(x=SerDep, y=DUSP6, label = rownames(pol_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor pol", subtitle = stats)

## TEs
### ATF4 (sanity check)
TEs_plots.cor <-rcorr(TEs_plots$SerDep, TEs_plots$ATF4, type = "spearman")
stats <- paste("Spearman: R = ", round(TEs_plots.cor$r["x","y"], 2),
               ", pval = ", round(TEs_plots.cor$P["x","y"], 4), sep = "")

ggplot(TEs_plots, aes(x=SerDep, y=ATF4, label = rownames(TEs_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor residuals", subtitle = stats)

### DUSP6 
TEs_plots.cor <-rcorr(TEs_plots$SerDep, TEs_plots$DUSP6, type = "spearman")
stats <- paste("Spearman: R = ", round(TEs_plots.cor$r["x","y"], 2),
               ", pval = ", round(TEs_plots.cor$P["x","y"], 4), sep = "")

ggplot(TEs_plots, aes(x=SerDep, y=DUSP6, label = rownames(TEs_plots))) + geom_point() + 
  geom_smooth(method=lm) + theme_bw() + geom_text_repel() + labs(title="Tumor residuals", subtitle = stats)


# Find IC3 indicator
## PHGDH & CBS
ggplot(gene_plots, aes(x=CBS, y=PHGDH, color= SerDep, label = rownames(gene_plots))) + 
  geom_point() + 
#  scale_x_continuous(trans='log10') +
#  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "PHGDH/CBS as IC3 indicator")

## Only CBS
ggplot(gene_plots, aes(x=CBS, y=SerDep, color= SerDep, label = rownames(gene_plots))) + 
  geom_point() + 
#  scale_x_continuous(trans='log10') +
#  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "CBS as IC3 indicator")


## Only PHGDH
ggplot(gene_plots, aes(x=PHGDH, y=SerDep, color= SerDep, label = rownames(gene_plots))) + 
  geom_point() + 
  #  scale_x_continuous(trans='log10') +
  #  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "PHGDH as IC3 indicator")


## Only IL18
ggplot(gene_plots, aes(x=IL18, y=SerDep, color= SerDep, label = rownames(gene_plots))) + 
  geom_point() + 
  #  scale_x_continuous(trans='log10') +
  #  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "IL18 as IC3 indicator")

## IL18 and PHGDH
ggplot(gene_plots, aes(x=IL18, y=PHGDH, color= SerDep, label = rownames(gene_plots))) + 
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

gene_plots.cor <-rcorr(gene_plots$SerDep, c(gene_plots$comIL18PHGDH1), type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")
print(stats)
ggplot(gene_plots, aes(x=comIL18PHGDH1, y=SerDep, color= SerDep, label = rownames(gene_plots))) + 
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

gene_plots.cor <-rcorr(gene_plots$SerDep, c(gene_plots$comIL18PHGDH2), type = "spearman")
stats <- paste("Spearman: R = ", round(gene_plots.cor$r["x","y"], 2),
               ", pval = ", round(gene_plots.cor$P["x","y"], 4), sep = "")
print(stats)
ggplot(gene_plots, aes(x=comIL18PHGDH2, y=SerDep, color= SerDep, label = rownames(gene_plots))) + 
  geom_point() + 
  #  scale_x_continuous(trans='log10') +
  #  scale_y_continuous(trans='log10') +
  geom_smooth(method=lm) +
  scale_colour_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  #geom_text_repel() +
  labs(title="28PDX Tumor", subtitle = "combined IL18 PHGDH mean(IL18, -PHGDH) as IC3 indicator")

